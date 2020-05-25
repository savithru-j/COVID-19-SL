#include <iostream>
#include <chrono>

#include "CovidSim.h"

int
main()
{
  ModelParams params(75, 14);

  Population pop0(21323733, 5, 1);

  auto t0 = std::chrono::high_resolution_clock::now();

  std::vector<Population> pop_hist;
  pop_hist = predictModel(params, pop0);

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;

  for (int i = 0; i < 10; i++)
    std::cout << uniformRand() << std::endl;

	ObservedPopulation pop_SL("csv_data/srilanka.txt");

//	Vector x = std::vector<double>{0.4, 5, 2.7, -0.3, 3.6};
//	Matrix A = getHaarMatrix(x.size());
//	Vector y = A*x;

	getOptimalParameters(pop_SL);

	return 0;
}


std::vector<Population>
predictModel(const ModelParams& params, const Population& pop_init)
{
  std::vector<Population> population_hist;

  const int nt = params.nt_hist + params.nt_pred;
  const int nt_sub = 1.0/params.dt;

  Population pop = pop_init;

  population_hist.reserve(nt);
  population_hist.push_back(pop);

  for (int t = 0; t < nt-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);

    pop.report(params, t); //report once daily

    population_hist.push_back(pop); //save solution daily
  }

  return population_hist;
}

ModelParams
getOptimalParameters(const ObservedPopulation& observed_pop)
{
  f_eval_count = 0;
  const int T_end_buffer = 5; //no. of non-optimized days at end

  const int nt = observed_pop.getNumDays();
  const int nt_opt = nt - T_end_buffer; //number of days to optimize parameters for

  const auto param_bounds = getParameterBounds(nt_opt);
  Vector param_vec(param_bounds.size());

  double jump_energy = 1.0;
  randomizeParameterVector(param_bounds, jump_energy, param_vec);

  ModelParams params(nt_opt, T_end_buffer);
  copyVector2Param(param_vec, params);

  Vector param_vec0 = param_vec;
  Vector param_vec_opt = param_vec;

  Population pop_init(observed_pop.N, 5, 1); //TODO: set initial E0 and Rd

  Matrix reg_matrix = getHaarMatrix(nt_opt);

  double cost0 = getCost(params, pop_init, observed_pop);
  double err_scaling = 1.0 / cost0;

  cost0 = getCost(params, pop_init, observed_pop, reg_matrix, err_scaling);
//  double cost = cost0;
//  Vector cost_grad = getCostGradient(params, pop_init, observed_pop, param_bounds, reg_matrix, err_scaling);

  double cost_rel_min = 1.0; //cost / cost0;
  const double cost_reduction_tol = 1e-4;
  const double min_eta = 1e-6;
  const int max_iter_per_pass = 10;
  const int max_passes = 5;

  std::cout << std::scientific << std::setprecision(4);

  for (int pass = 0; pass < max_passes; ++pass)
  {
    std::cout << std::endl << "Pass " << pass << " ----------------" << std::endl;

    double cost = getCost(params, pop_init, observed_pop, reg_matrix, err_scaling);

    double cost_prev = cost;
    int stalled_iter = 0;

    int iter = 0;
    for (iter = 0; iter < max_iter_per_pass; ++iter)
    {
      Vector cost_grad = getCostGradient(params, pop_init, observed_pop, param_bounds, reg_matrix, err_scaling);

      double eta0 = limitUpdate(param_bounds, param_vec, cost_grad);
      std::cout << "  Iter: " << iter << ", cost: " << (cost/cost0) << ", eta0: " << eta0 << std::endl;

      if (eta0 < 1e-13)
        break;

      for (int i = 0; i < param_vec.m(); ++i)
      {
        param_vec0[i] = param_vec[i];
        cost_grad[i] *= eta0;
      }

      double eta = 1.0;
      while (eta >= min_eta)
      {
        //Update param_vec based on (scaled) update vector
        for (int i = 0; i < param_vec.m(); ++i)
          param_vec[i] = param_vec0[i] - eta*cost_grad[i];

        copyVector2Param(param_vec, params);

        //Evaluate new cost
        double cost_new = getCost(params, pop_init, observed_pop, reg_matrix, err_scaling);
//        std::cout << "    eta: " << eta << ", cost: " << (cost_new/cost0) << std::endl;

        if (cost_new < cost)
        {
          cost = cost_new;
          break;
        }

        eta /= 2.0;
      } //linesearch

      //Check if residual has stalled
      if ((cost_prev - cost) < cost_prev*cost_reduction_tol)
        stalled_iter++;
      else {
        cost_prev = cost;
        stalled_iter = 0;
      }

      if (stalled_iter == 5) //Abort if residual has stalled for too many iterations.
      {
        std::cout << "Descent stalled." << std::endl;
        break;
      }

    } //gradient descent

    double cost_rel = cost/cost0;
    std::cout << "Num_iter: " << iter << ", cost: " << (cost/cost0) << std::endl;

    if (cost_rel < cost_rel_min)
    {
      cost_rel_min = cost_rel;
      param_vec_opt = param_vec;
    }

    jump_energy = 1.0 - pass / (double) (max_passes - 1.0);
    std::cout << "Jump energy: " << jump_energy << std::endl;

    param_vec = param_vec_opt;
    randomizeParameterVector(param_bounds, jump_energy, param_vec);
    copyVector2Param(param_vec, params);

  } //passes

  std::cout << "Minimum cost: " << cost_rel_min << std::endl;
  std::cout << "Evaluation count: " << f_eval_count << std::endl;

  return params;
}

double getCost(const ModelParams& params, const Population& pop_init,
               const ObservedPopulation& observed_pop, const Matrix& reg_matrix,
               const double& scaling)
{
  auto pop_hist = predictModel(params, pop_init);

  const double wC = 1.0, wR = 1.0, wF = 1.0;
  const int nt = params.nt_hist + params.nt_pred;

  double cost = 0.0;
  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - observed_pop.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - observed_pop.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - observed_pop.deaths[i]);

    cost += wC*err_total*err_total + wR*err_recov*err_recov + wF*err_fatal*err_fatal;
  }

  cost *= scaling;

  //Add regularization terms
  if (reg_matrix.m() > 0 && reg_matrix.n() > 0)
  {
    double reg_betaN = 0.0;
    double reg_ce = 0.0;
    double reg_c1 = 0.0;

    for (int i = 0; i < reg_matrix.m(); ++i)
    {
      double coeff_betaN = 0.0;
      double coeff_ce = 0.0;
      double coeff_c1 = 0.0;

      for (int j = 0; j < reg_matrix.n(); ++j)
      {
        if (j < params.betaN.m())
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN(j);
          coeff_ce += reg_matrix(i,j) * params.ce(j);
          coeff_c1 += reg_matrix(i,j) * params.c1(j);
        }
        else
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN.back();
          coeff_ce += reg_matrix(i,j) * params.ce.back();
          coeff_c1 += reg_matrix(i,j) * params.c1.back();
        }
      }

      //Take L1-norm of coefficient vector
      reg_betaN += std::abs(coeff_betaN);
      reg_ce += std::abs(coeff_ce);
      reg_c1 += std::abs(coeff_c1);
    }

    const double wgt = 1e-1;
    const double reg_cost = wgt*reg_betaN + wgt*reg_ce + wgt*reg_c1;
    cost += reg_cost;

//    if (std::isnan(cost))
//      throwError("Found NaN!");
  }

//  if (isNaN(cost))
//  {
//    console.log("getOptimizationCost: found NaN");
//    return {};
//  }
  f_eval_count++;
  return cost;
}

Vector getCostGradient(const ModelParams& params_orig, const Population& pop_init,
                       const ObservedPopulation& observed_pop, const std::vector<ParamBound>& bounds,
                       const Matrix& reg_matrix, const double& scaling)
{
  Vector param_vec(bounds.size());
  Vector grad(bounds.size());

  ModelParams params = params_orig; //create a copy
  copyParam2Vector(params, param_vec);

  for (int j = 0; j < param_vec.m(); ++j)
  {
    double delta = bounds[j].step;

    //Compute finite difference
    param_vec[j] += delta;
    copyVector2Param(param_vec, params);
    double fp = getCost(params, pop_init, observed_pop, reg_matrix, scaling);

    param_vec[j] -= 2*delta;
    copyVector2Param(param_vec, params);
    double fm = getCost(params, pop_init, observed_pop, reg_matrix, scaling);

    param_vec[j] += delta;

    grad[j] = (fp - fm)/(2*delta);
  }
  return grad;
}

void copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int nt = params.nt_hist;
  const int m = 3*nt + 10;
  if (v.m() != m)
    v.resize(m);

  for (int i = 0; i < nt; ++i)
  {
    v[       i] = params.betaN[i];
    v[  nt + i] = params.ce[i];
    v[2*nt + i] = params.c1[i];
  }
  const int off = 3*nt;
  v[off  ] = params.T_incub0;
  v[off+1] = params.T_incub1;
  v[off+2] = params.T_asympt;
  v[off+3] = params.T_mild;
  v[off+4] = params.T_severe;
  v[off+5] = params.T_icu;
  v[off+6] = params.f;
  v[off+7] = params.frac_recover_I1;
  v[off+8] = params.frac_recover_I2;
  v[off+9] = params.CFR;
}

void copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 3*nt + 10;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

//  for (int i = 0; i < v.m(); ++i)
//    if (std::isnan(v[i]))
//      throwError("Found NaN!");

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[       i];
    params.ce[i]    = v[  nt + i];
    params.c1[i]    = v[2*nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = v[  nt-1];
    params.ce[i]    = v[2*nt-1];
    params.c1[i]    = v[3*nt-1];
  }

  const int off = 3*nt;
  params.T_incub0        = v[off  ];
  params.T_incub1        = v[off+1];
  params.T_asympt        = v[off+2];
  params.T_mild          = v[off+3];
  params.T_severe        = v[off+4];
  params.T_icu           = v[off+5];
  params.f               = v[off+6];
  params.frac_recover_I1 = v[off+7];
  params.frac_recover_I2 = v[off+8];
  params.CFR             = v[off+9];
}

std::vector<ParamBound>
getParameterBounds(int nt)
{
  const int m = 3*nt + 10;
  std::vector<ParamBound> bounds(m);

  for (int i = 0; i < nt; ++i)
  {
    bounds[       i] = ParamBound(0, 1, 1e-4); //betaN
    bounds[  nt + i] = ParamBound(0, 1, 1e-4); //ce
    bounds[2*nt + i] = ParamBound(0, 1, 1e-4); //c1
  }
  const int off = 3*nt;
  bounds[off  ] = ParamBound(1, 10, 1e-4); //T_incub0
  bounds[off+1] = ParamBound(1, 10, 1e-4); //T_incub1
  bounds[off+2] = ParamBound(1, 10, 1e-4); //T_asympt
  bounds[off+3] = ParamBound(1, 10, 1e-4); //T_mild
  bounds[off+4] = ParamBound(1, 10, 1e-4); //T_severe
  bounds[off+5] = ParamBound(1, 10, 1e-4); //T_icu
  bounds[off+6] = ParamBound(0, 1, 1e-4); //f
  bounds[off+7] = ParamBound(0, 1, 1e-4); //frac_recover_I1
  bounds[off+8] = ParamBound(0, 1, 1e-4); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 1e-5); //CFR

  return bounds;
}

double uniformRand(double min, double max)
{
//  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(2.0); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> uniform_rand(min, max);
  return uniform_rand(gen);
}

void
randomizeParameterVector(const std::vector<ParamBound>& bounds, const double energy, Vector& param_vec)
{
  for (int i = 0; i < param_vec.m(); ++i)
  {
    const double delta = energy * (bounds[i].max - bounds[i].min) * uniformRand(-1, 1);
    param_vec[i] = std::min(std::max(param_vec[i] + delta, bounds[i].min), bounds[i].max);
  }
}

double
limitUpdate(const std::vector<ParamBound>& bounds, const Vector& param_vec, Vector& dparam_vec)
{
  double eta = 1.0;

  for (int i = 0; i < param_vec.m(); ++i)
  {
    if ( (param_vec[i] == bounds[i].min && -dparam_vec[i] < 0.0) ||
         (param_vec[i] == bounds[i].max && -dparam_vec[i] > 0.0) )
    {
       dparam_vec[i] = 0.0;
    }
    else
    {
      double new_val = param_vec[i] - dparam_vec[i];
      if (new_val < bounds[i].min)
        eta = std::min(eta, (param_vec[i] - bounds[i].min)/dparam_vec[i]);
      else if (new_val > bounds[i].max)
        eta = std::min(eta, (param_vec[i] - bounds[i].max)/dparam_vec[i]);
    }
  }
  return eta;
}

Matrix getHaarMatrix(int m)
{
  const int num_levels = std::ceil(std::log2(m));
  const int n = std::pow(2, num_levels);
//  if (std::pow(2,num_levels) != n)
//    throw("N is not a power of 2!");

  Matrix A(n, n, 0.0);
  const double inv_sqrtn = 1.0 / std::sqrt(n);

  for (int j = 0; j < n; ++j)
    A(0, j) = inv_sqrtn;

  for (int k = 1; k < n; ++k)
  {
    int p = std::floor(std::log(k) / std::log(2.0));
    double k1 = std::pow(2.0, p);
    double k2 = k1 * 2.0;
    double q = k - k1;
    double t1 = n/k1;
    double t2 = n/k2;
    double tmp = std::pow(2.0, p/2.0) * inv_sqrtn;
    for (int i = 0; i < t2; ++i)
    {
      A(k, q*t1 + i) = tmp;
      A(k, + q*t1 + i + t2) = -tmp;
    }
  }

  return A;
}
