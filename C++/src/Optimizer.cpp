#include <iostream>
#include <iomanip>
#include <random>

#include <nlopt.hpp>

#include "Optimization.h"
#include "Simulator.h"

Optimizer::Optimizer(const ObservedPopulation& pop_observed_,
                     const Population& pop_init_, double reg_weight_) :
    pop_observed(pop_observed_), pop_init(pop_init_), reg_weight(reg_weight_),
    nt_opt(pop_observed.getNumDays() - t_buffer),
    params(nt_opt, t_buffer)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  if (reg_weight != 0.0)
    reg_matrix = getHaarMatrix(nt_opt);

  param_bounds = getParameterBounds(nt_opt);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}

void
Optimizer::optimizeParameters()
{
  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;

  f_eval_count = 0;

  double jump_energy = 1.0;

  Vector param_vec0 = param_vec;
  optimal_param_vec.fill(param_vec);

  double cost0 = -1;
  cost_rel_min.fill(1.0);

  Vector cost_grad(nDim(), 0.0); //storage vector for cost gradient

  std::cout << std::scientific << std::setprecision(4);

  for (int pass = 0; pass < max_passes; ++pass)
  {
    std::cout << std::endl << "Pass " << pass << " ----------------" << std::endl;

    //max_iter_per_pass = 50 + (pass/(double)(max_passes - 1)) * (500 - 50);

    auto cost = getCost();
    if (cost0 == -1)
      cost0 = cost.first; //save initial cost

    double cost_prev = cost.first;
    int stalled_iter = 0;

    int iter = 0;
    for (iter = 0; iter < max_iter_per_pass; ++iter)
    {
      getCostGradient(cost_grad);

      double eta0 = limitUpdate(cost_grad);
      std::cout << "  Iter: " << iter << ", cost: " << (cost.first/cost0) << ", eta0: " << eta0 << std::endl;

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
        auto cost_new = getCost();
//        std::cout << "    eta: " << eta << ", cost: " << (cost_new/cost0) << std::endl;

        if (cost_new.first < cost.first)
        {
          cost = cost_new;
          break;
        }

        eta /= 2.0;
      } //linesearch

      //Check if residual has stalled
      if ((cost_prev - cost.first) < cost_prev*cost_reduction_tol)
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

    if (cost_rel < cost_rel_min.back())
    {
      updateOptimalSolution(cost_rel, cost.second, param_vec);
    }

    jump_energy = 1.0 - pass / (double) (max_passes - 1.0);
    std::cout << "Jump energy: " << jump_energy << std::endl;

    param_vec = optimal_param_vec[0]; //goto current optimal solution
    randomizeParameters(jump_energy);

  } //passes

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_rel_min << std::endl;
  std::cout << "Minimum sub-costs [sol-error, beta-reg, c0-reg, c1-reg]:" << std::endl;
  for (std::size_t i = 0; i < NUM_RESULTS; i++)
      std::cout << "  " << sub_costs_min[i] << std::endl;

  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;

//  for (std::size_t i = 0; i < NUM_RESULTS; i++)
//    std::cout << "Optimal params[" << i << "]: " << param_vec_opt[i] << std::endl << std::endl;
}

void
Optimizer::optimizeParametersNLOPT()
{
  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;

  f_eval_count = 0;

//  const auto param_bounds = getParameterBounds(nt_opt);
//  std::size_t ndim = param_bounds.size();
//  Vector param_vec(ndim);

  //nlopt optimizer object
  nlopt::opt opt(nlopt::LD_LBFGS, nDim());

  opt.set_min_objective( getCostNLOPT, reinterpret_cast<void*>(this) );

  //stop when every parameter changes by less than the tolerance multiplied by the absolute value of the parameter.
  opt.set_xtol_rel(1e-6);

  //stop when the objective function value changes by less than the tolerance multiplied by the absolute value of the function value
  opt.set_ftol_rel(1e-8);

  //stop when the maximum number of function evaluations is reached
  opt.set_maxeval(100);

//  opt.set_lower_bounds({2, -1000});

  std::vector<double> x = param_vec.getDataVector();
  double f_opt;

  try
  {
    opt.optimize(x, f_opt);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    throwError("NLopt failed!");
  }

  for (int i = 0; i < NUM_RESULTS; ++i)
    optimal_param_vec[i] = x;

  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
  std::cout << "Optimal cost: " << f_opt << std::endl;
}

double
Optimizer::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  Optimizer* opt = reinterpret_cast<Optimizer*>(data);

  opt->param_vec = x;
  copyVector2Param(opt->param_vec, opt->params);

  double cost = opt->getCost().first;

  if (!grad.empty())
    opt->getCostGradient(grad);

  std::cout << "cost: " << cost << std::endl;

  return cost;
//  if (!grad.empty())
//  {
//    grad[0] = 2*(x[0]);
//    grad[1] = 2*(x[1]);
//  }
//  double f = x[0]*x[0] + x[1]*x[1];
//  std::cout << x[0] << ", " << x[1] << ": " << f << std::endl;
//  return f;
}


std::pair<double, std::array<double, 4>>
Optimizer::getCost()
{
  auto pop_hist = predictModel(params, pop_init);

  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,4> sub_costs = {0.0, 0.0, 0.0, 0.0};

  double err_sq_total = 0.0, err_sq_recov = 0.0, err_sq_fatal = 0.0;
  double sq_total = 0.0, sq_recov = 0.0, sq_fatal = 0.0;

  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - pop_observed.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - pop_observed.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - pop_observed.deaths[i]);

    err_sq_total += err_total*err_total;
    err_sq_recov += err_recov*err_recov;
    err_sq_fatal += err_fatal*err_fatal;

    sq_total += pop_observed.confirmed[i]*pop_observed.confirmed[i];
    sq_recov += pop_observed.recovered[i]*pop_observed.recovered[i];
    sq_fatal += pop_observed.deaths[i]*pop_observed.deaths[i];
  }

  const double wC = 1.0, wR = 0.0, wF = 1.0;
  sub_costs[0] = wC*(err_sq_total/sq_total) + wR*(err_sq_recov/sq_recov) + wF*(err_sq_fatal/sq_fatal);

  //Add regularization terms
  if (reg_weight != 0.0)
  {
    for (int i = 0; i < reg_matrix.m(); ++i)
    {
      double coeff_betaN = 0.0;
      double coeff_c0 = 0.0;
      double coeff_c1 = 0.0;

      for (int j = 0; j < reg_matrix.n(); ++j)
      {
        if (j < params.betaN.m())
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN(j);
          coeff_c0 += reg_matrix(i,j) * params.c0(j);
          coeff_c1 += reg_matrix(i,j) * params.c1(j);
        }
        else
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN.back();
          coeff_c0 += reg_matrix(i,j) * params.c0.back();
          coeff_c1 += reg_matrix(i,j) * params.c1.back();
        }
      }

      //Take L1-norms of coefficient vectors
      sub_costs[1] += std::abs(coeff_betaN);
      sub_costs[2] += std::abs(coeff_c0);
      sub_costs[3] += std::abs(coeff_c1);
    }
    sub_costs[1] *= reg_weight;
    sub_costs[2] *= reg_weight;
    sub_costs[3] *= reg_weight;
  } //if-regularize

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2] + sub_costs[3];

  f_eval_count++;
  return {cost, sub_costs};
}

void
Optimizer::getCostGradient(std::vector<double>& grad)
{
  ModelParams params_orig = params; //create a copy of the current parameters

  for (int j = 0; j < param_vec.m(); ++j)
  {
    const double delta = param_bounds[j].step;

    //Compute finite difference
    param_vec[j] += delta;
    copyVector2Param(param_vec, params);
    double fp = getCost().first;

    param_vec[j] -= 2*delta;
    copyVector2Param(param_vec, params);
    double fm = getCost().first;

    param_vec[j] += delta;

    grad[j] = (fp - fm)/(2*delta);
  }

  params = params_orig; //Restore current parameters
  copyParam2Vector(params, param_vec);
}

void
Optimizer::updateOptimalSolution(
    const double& cost_rel, const std::array<double,4>& sub_costs, const Vector& param_vec_cur)
{
  int min_ind = 0;
  for (min_ind = 0; min_ind < NUM_RESULTS; ++min_ind)
    if (cost_rel < cost_rel_min[min_ind])
      break;

  for (int i = NUM_RESULTS-1; i > min_ind; i--)
  {
    cost_rel_min[i] = cost_rel_min[i-1];
    sub_costs_min[i] = sub_costs_min[i-1];
    optimal_param_vec[i] = optimal_param_vec[i-1];
  }
  cost_rel_min[min_ind] = cost_rel;
  sub_costs_min[min_ind] = sub_costs;
  optimal_param_vec[min_ind] = param_vec_cur;
}

double
Optimizer::limitUpdate(Vector& dparam_vec)
{
  double eta = 1.0;

  for (int i = 0; i < param_vec.m(); ++i)
  {
    if ( (param_vec[i] == param_bounds[i].min && -dparam_vec[i] < 0.0) ||
         (param_vec[i] == param_bounds[i].max && -dparam_vec[i] > 0.0) )
    {
       dparam_vec[i] = 0.0;
    }
    else
    {
      double new_val = param_vec[i] - dparam_vec[i];
      if (new_val < param_bounds[i].min)
        eta = std::min(eta, (param_vec[i] - param_bounds[i].min)/dparam_vec[i]);
      else if (new_val > param_bounds[i].max)
        eta = std::min(eta, (param_vec[i] - param_bounds[i].max)/dparam_vec[i]);
    }
  }
  return eta;
}


void copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int nt = params.nt_hist;
  const int m = 3*nt + 10;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");
//  if (v.m() != m)
//    v.resize(m);

  for (int i = 0; i < nt; ++i)
  {
    v[       i] = params.betaN[i];
    v[  nt + i] = params.c0[i];
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

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[       i];
    params.c0[i]    = v[  nt + i];
    params.c1[i]    = v[2*nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = v[  nt-1];
    params.c0[i]    = v[2*nt-1];
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

  const double delta = 1e-4;

  for (int i = 0; i < nt; ++i)
  {
    bounds[       i] = ParamBound(0, 2, delta); //betaN
    bounds[  nt + i] = ParamBound(0, 1, delta); //c0
    bounds[2*nt + i] = ParamBound(0, 1, delta); //c1
  }
  const int off = 3*nt;
  bounds[off  ] = ParamBound(2.9, 3.1, delta); //T_incub0
  bounds[off+1] = ParamBound(1.9, 2.1, delta); //T_incub1
  bounds[off+2] = ParamBound(5.9, 6.1, delta); //T_asympt
  bounds[off+3] = ParamBound(5.9, 6.1, delta); //T_mild
  bounds[off+4] = ParamBound(3.9, 4.1, delta); //T_severe
  bounds[off+5] = ParamBound(9.9, 10.1, delta); //T_icu
  bounds[off+6] = ParamBound(0.29, 0.31, delta); //f
  bounds[off+7] = ParamBound(0.5, 0.9, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.5, 0.9, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR

  return bounds;
}

double uniformRand(double min, double max)
{
//  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(2.0); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> uniform_rand(min, max);
  return uniform_rand(gen);
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
