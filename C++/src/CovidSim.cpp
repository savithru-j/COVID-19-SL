#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>

#include "CovidSim.h"

int
main()
{
  std::string country = "srilanka";

	ObservedPopulation pop_observed("csv_data/" + country + ".txt");

	Population pop_init(pop_observed.N, 5, 1);

	std::string folder_path = "results";
	mkdir(folder_path.c_str(), 0777);

	std::string filepath_opt_params = folder_path + "/" + country + "_params.txt";
  std::ofstream file_opt_params(filepath_opt_params);
  if (!file_opt_params.good())
    throwError("Cannot open file to write - " + filepath_opt_params);

  std::string filepath_predictions = folder_path + "/" + country + "_prediction.txt";
  std::ofstream file_predictions(filepath_predictions);
  if (!file_predictions.good())
    throwError("Cannot open file to write - " + filepath_predictions);

  auto t0 = std::chrono::high_resolution_clock::now();

  OptimizationInfo optinfo = optimizeParameters(pop_observed, pop_init);

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;


  file_opt_params << std::scientific << std::setprecision(6);
  for (int i = 0; i < optinfo.param_vec_opt[0].m(); ++i)
  {
    for (std::size_t j = 0; j < optinfo.param_vec_opt.size()-1; ++j)
      file_opt_params << optinfo.param_vec_opt[j][i] << ", ";
    file_opt_params << optinfo.param_vec_opt.back()[i] << std::endl;
  }

  const int nt_hist = pop_observed.getNumDays();

  std::array<std::vector<Population>, OptimizationInfo::NUM_RESULTS> predictions;
  for (std::size_t i = 0; i < predictions.size(); ++i)
  {
    ModelParams params(nt_hist - 5, 5+14);
    copyVector2Param(optinfo.param_vec_opt[i], params);
    predictions[i] = predictModel(params, pop_init);
  }

  file_predictions << std::scientific << std::setprecision(6);
  for (std::size_t i = 0; i < predictions[0].size(); ++i)
  {
    for (std::size_t j = 0; j < predictions.size(); ++j)
    file_predictions << predictions[j][i].getNumReported() << ", "
                     << predictions[j][i].getNumRecoveredReported() << ", "
                     << predictions[j][i].getNumFatalReported() << ", ";
    file_predictions << std::endl;
  }

  file_opt_params.close();
  file_predictions.close();

//  ModelParams params(nt_hist - 5, 0);
//  copyVector2Param(optinfo.param_vec_opt[0], params);
//  Vector x(128);
//  for (int i = 0; i < 128; ++i)
//  {
//    if (i < params.betaN.m())
//      x[i] = params.betaN[i];
//    else
//      x[i] = params.betaN.back();
//  }
//
//  Vector coeff = optinfo.reg_matrix*x;
//
//  std::cout << "coeff: " << coeff << std::endl;

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

OptimizationInfo
optimizeParameters(const ObservedPopulation& pop_observed, const Population& pop_init)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  OptimizationInfo optinfo;

  const int nt = pop_observed.getNumDays();
  const int nt_opt = nt - OptimizationInfo::T_END_BUFFER; //number of days to optimize parameters for

  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;

  optinfo.reg_matrix = getHaarMatrix(nt_opt);

  const auto param_bounds = getParameterBounds(nt_opt);
  Vector param_vec(param_bounds.size());

  double jump_energy = 1.0;
  randomizeParameterVector(param_bounds, jump_energy, param_vec);

  ModelParams params(nt_opt, OptimizationInfo::T_END_BUFFER);
  copyVector2Param(param_vec, params);

  Vector param_vec0 = param_vec;
  optinfo.param_vec_opt.fill(param_vec);

  optinfo.regularize = false; //turn off regularization to get initial raw cost
  double cost0 = getCost(params, pop_init, pop_observed, optinfo).first; //raw cost without regularization or scaling

  optinfo.regularize = true;
  optinfo.cost_error_scaling = 1.0 / cost0;
  cost0 = getCost(params, pop_init, pop_observed, optinfo).first;

  optinfo.cost_rel_min.fill(1.0);
  optinfo.cost_raw_min.fill(1.0);

  std::cout << std::scientific << std::setprecision(4);

  for (int pass = 0; pass < OptimizationInfo::MAX_PASSES; ++pass)
  {
    std::cout << std::endl << "Pass " << pass << " ----------------" << std::endl;

    auto cost = getCost(params, pop_init, pop_observed, optinfo);

    double cost_prev = cost.first;
    int stalled_iter = 0;

    int iter = 0;
    for (iter = 0; iter < OptimizationInfo::MAX_ITER_PER_PASS; ++iter)
    {
      Vector cost_grad = getCostGradient(params, pop_init, pop_observed, param_bounds, optinfo);

      double eta0 = limitUpdate(param_bounds, param_vec, cost_grad);
      std::cout << "  Iter: " << iter << ", cost: " << (cost.first/cost0) << ", eta0: " << eta0 << std::endl;

      if (eta0 < 1e-13)
        break;

      for (int i = 0; i < param_vec.m(); ++i)
      {
        param_vec0[i] = param_vec[i];
        cost_grad[i] *= eta0;
      }

      double eta = 1.0;
      while (eta >= OptimizationInfo::MIN_ETA)
      {
        //Update param_vec based on (scaled) update vector
        for (int i = 0; i < param_vec.m(); ++i)
          param_vec[i] = param_vec0[i] - eta*cost_grad[i];

        copyVector2Param(param_vec, params);

        //Evaluate new cost
        auto cost_new = getCost(params, pop_init, pop_observed, optinfo);
//        std::cout << "    eta: " << eta << ", cost: " << (cost_new/cost0) << std::endl;

        if (cost_new.first < cost.first)
        {
          cost = cost_new;
          break;
        }

        eta /= 2.0;
      } //linesearch

      //Check if residual has stalled
      if ((cost_prev - cost.first) < cost_prev*OptimizationInfo::COST_REDUCTION_TOL)
        stalled_iter++;
      else {
        cost_prev = cost.first;
        stalled_iter = 0;
      }

      if (stalled_iter == 5) //Abort if residual has stalled for too many iterations.
      {
        std::cout << "Descent stalled." << std::endl;
        break;
      }

    } //gradient descent

    double cost_rel = cost.first/cost0;
    std::cout << "Num_iter: " << iter << ", cost: " << cost_rel << std::endl;

    if (cost_rel < optinfo.cost_rel_min.back())
    {
      bool tmp_regularize = optinfo.regularize;
      double tmp_scaling = optinfo.cost_error_scaling;

      optinfo.regularize = false;
      optinfo.cost_error_scaling = 1.0;
      double cost_raw = getCost(params, pop_init, pop_observed, optinfo).first; //raw cost without regularization or scaling

      //revert
      optinfo.regularize = tmp_regularize;
      optinfo.cost_error_scaling = tmp_scaling;

      optinfo.updateOptimalSolution(cost_rel, cost_raw, cost.second, param_vec);
    }

    jump_energy = 1.0 - pass / (double) (OptimizationInfo::MAX_PASSES - 1.0);
    std::cout << "Jump energy: " << jump_energy << std::endl;

    param_vec = optinfo.param_vec_opt[0];
    randomizeParameterVector(param_bounds, jump_energy, param_vec);
    copyVector2Param(param_vec, params);

  } //passes

  std::cout << std::endl;
  std::cout << "Error scaling: " << optinfo.cost_error_scaling << std::endl;

  std::cout << "Minimum costs (relative): " << optinfo.cost_rel_min << std::endl;
  std::cout << "Minimum costs (raw): " << optinfo.cost_raw_min << std::endl;
  std::cout << "Minimum sub-costs [sol-error, beta-reg, ce-reg, c1-reg]:" << std::endl;
  for (std::size_t i = 0; i < OptimizationInfo::NUM_RESULTS; i++)
      std::cout << "  " << optinfo.sub_costs_min[i] << std::endl;

  std::cout << "Cost function evaluation count: " << optinfo.f_eval_count << std::endl;

//  for (std::size_t i = 0; i < NUM_RESULTS; i++)
//    std::cout << "Optimal params[" << i << "]: " << param_vec_opt[i] << std::endl << std::endl;

  return optinfo;
}

std::pair<double, std::array<double, 4>>
getCost(const ModelParams& params, const Population& pop_init,
        const ObservedPopulation& observed_pop, OptimizationInfo& optinfo)
{
  auto pop_hist = predictModel(params, pop_init);

  const double wC = 1.0, wR = 1.0, wF = 10.0;
  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,4> sub_costs = {0.0, 0.0, 0.0, 0.0};

  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - observed_pop.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - observed_pop.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - observed_pop.deaths[i]);

    sub_costs[0] += wC*err_total*err_total + wR*err_recov*err_recov + wF*err_fatal*err_fatal;
  }

  sub_costs[0] *= optinfo.cost_error_scaling;

  //Add regularization terms
  if (optinfo.regularize)
  {
    for (int i = 0; i < optinfo.reg_matrix.m(); ++i)
    {
      double coeff_betaN = 0.0;
      double coeff_ce = 0.0;
      double coeff_c1 = 0.0;

      for (int j = 0; j < optinfo.reg_matrix.n(); ++j)
      {
        if (j < params.betaN.m())
        {
          coeff_betaN += optinfo.reg_matrix(i,j) * params.betaN(j);
          coeff_ce += optinfo.reg_matrix(i,j) * params.ce(j);
          coeff_c1 += optinfo.reg_matrix(i,j) * params.c1(j);
        }
        else
        {
          coeff_betaN += optinfo.reg_matrix(i,j) * params.betaN.back();
          coeff_ce += optinfo.reg_matrix(i,j) * params.ce.back();
          coeff_c1 += optinfo.reg_matrix(i,j) * params.c1.back();
        }
      }

      //Take L1-norms of coefficient vectors
      sub_costs[1] += std::abs(coeff_betaN);
      sub_costs[2] += std::abs(coeff_ce);
      sub_costs[3] += std::abs(coeff_c1);
    }
    const double wReg = 0.01;
    sub_costs[1] *= wReg;
    sub_costs[2] *= wReg;
    sub_costs[3] *= wReg;
  } //if-regularize

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2] + sub_costs[3];

  optinfo.f_eval_count++;
  return {cost, sub_costs};
}

Vector getCostGradient(const ModelParams& params_orig, const Population& pop_init,
                       const ObservedPopulation& pop_observed, const std::vector<ParamBound>& bounds,
                       OptimizationInfo& optinfo)
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
    double fp = getCost(params, pop_init, pop_observed, optinfo).first;

    param_vec[j] -= 2*delta;
    copyVector2Param(param_vec, params);
    double fm = getCost(params, pop_init, pop_observed, optinfo).first;

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
