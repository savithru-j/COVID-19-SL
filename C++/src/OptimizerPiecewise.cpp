#include <iostream>
#include <iomanip>

#include "OptimizerPiecewise.h"
#include "Simulator.h"

OptimizerPiecewise::OptimizerPiecewise(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                                 const Vector& quarantine_input, int interval_size_, bool linear_basis_,
                                 double wconf, double wrecov, double wfatal,
                                 int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), interval_size(interval_size_), linear_basis(linear_basis_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0, quarantine_input),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  param_bounds = getParameterBounds(nt_opt, interval_size, linear_basis);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}


void
OptimizerPiecewise::optimizeParametersNLOPT()
{
  std::cout << "Optimizing parameters for " << params.nt_hist << " days..." << std::endl;
  std::cout << "No. of parameters: " << nDim() << std::endl;

  f_eval_count = 0;

  //nlopt optimizer object
  nlopt::opt opt(nlopt::LD_LBFGS, nDim());

  opt.set_min_objective( getCostNLOPT, reinterpret_cast<void*>(this) );

  //stop when every parameter changes by less than the tolerance multiplied by the absolute value of the parameter.
  opt.set_xtol_rel(1e-6);

  //stop when the objective function value changes by less than the tolerance multiplied by the absolute value of the function value
  opt.set_ftol_rel(1e-8);

  //stop when the maximum number of function evaluations is reached
  opt.set_maxeval(max_iter_per_pass);

  std::vector<double> lower_bounds(nDim()), upper_bounds(nDim());
  for (std::size_t i = 0; i < param_bounds.size(); ++i)
  {
    lower_bounds[i] = param_bounds[i].min;
    upper_bounds[i] = param_bounds[i].max;
  }
  opt.set_lower_bounds(lower_bounds);
  opt.set_upper_bounds(upper_bounds);

  cost_min.clear();
  sub_costs_min.clear();
  optimal_param_vec.clear();
  std::array<double,3> subcosts_dummy;

  std::cout << std::scientific << std::setprecision(4);

  for (int pass = 0; pass < max_passes; ++pass)
  {
    std::cout << std::endl << "Pass " << pass << "/" << max_passes << " ----------------" << std::endl;

    std::vector<double> x = param_vec.getDataVector();
    double cost_opt;
    nlopt::result result = nlopt::result::FAILURE;
    nlopt_iter = 0;

    try
    {
      result = opt.optimize(x, cost_opt);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      //throwError("NLopt failed!");
    }
    std::cout << "NLOPT result: " << getNLOPTResultDescription(result) << std::endl;
    std::cout << "Optimal cost: " << cost_opt << std::endl;

    updateOptimalSolution(cost_opt, subcosts_dummy, x);

    randomizeParameters();
    x = param_vec.getDataVector();
  }

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerPiecewise::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerPiecewise* opt = reinterpret_cast<OptimizerPiecewise*>(data);

  opt->param_vec = x;
  opt->copyVector2Param(opt->param_vec, opt->params);

  auto cost = opt->getCost();

  double grad_norm = 0.0;
  if (!grad.empty())
    grad_norm = opt->getCostGradient(grad);

  std::cout << "  Iter: " << opt->nlopt_iter << ", Cost: " << cost.first
            << " = [" << cost.second << "], Grad-norm: " << grad_norm << std::endl;
  opt->nlopt_iter++;

  return cost.first;
}


std::pair<double, std::array<double, 3>>
OptimizerPiecewise::getCost()
{
  auto pop_hist = predictModel(params, pop_init);

  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,3> sub_costs = {0.0, 0.0, 0.0};

  double err_sq_total = 0.0, err_sq_recov = 0.0, err_sq_fatal = 0.0;

#if 0 //L2
  double sq_total = 0.0, sq_recov = 0.0, sq_fatal = 0.0;
  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - pop_observed.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - pop_observed.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - pop_observed.deaths[i]);

    err_sq_total += err_total*err_total;
    err_sq_recov += err_recov*err_recov;
    err_sq_fatal += err_fatal*err_fatal;

    sq_total += (double)pop_observed.confirmed[i]*(double)pop_observed.confirmed[i];
    sq_recov += (double)pop_observed.recovered[i]*(double)pop_observed.recovered[i];
    sq_fatal += (double)pop_observed.deaths[i]   *(double)pop_observed.deaths[i];
  }
  sub_costs[0] = weight_conf*(err_sq_total/sq_total);
  sub_costs[1] = weight_recov*(err_sq_recov/sq_recov);
  sub_costs[2] = weight_fatal*(err_sq_fatal/sq_fatal);
#elif 0 //L2 of log
  for (int i = 1; i < nt; ++i)
  {
    double logA_obs = std::log(std::max(pop_observed.getNumActive(i), 1));
    double logR_obs = std::log(std::max(pop_observed.recovered[i], 1));
    double logF_obs = std::log(std::max(pop_observed.deaths[i], 1));

    double logA = std::log(std::max(pop_hist[i].getNumActiveReported(), 1.0));
    double logR = std::log(std::max(pop_hist[i].getNumRecoveredReported(), 1.0));
    double logF = std::log(std::max(pop_hist[i].getNumFatalReported(), 1.0));

    err_sq_total += (logA - logA_obs)*(logA - logA_obs);
    err_sq_recov += (logR - logR_obs)*(logR - logR_obs);
    err_sq_fatal += (logF - logF_obs)*(logF - logF_obs);
  }
  const double denom = (weight_conf + weight_recov + weight_fatal) * nt;
  sub_costs[0] = weight_conf * err_sq_total / denom;
  sub_costs[1] = weight_recov* err_sq_recov / denom;
  sub_costs[2] = weight_fatal* err_sq_fatal / denom;
#elif 1 //MLE
  constexpr double scaling = 1e-4;
  const double eps = 1e-5;
  for (int i = 1; i < nt; ++i)
  {
    double dC_obs = pop_observed.confirmed[i] - pop_observed.confirmed[i-1];
    double dR_obs = pop_observed.recovered[i] - pop_observed.recovered[i-1];
    double dF_obs = pop_observed.deaths[i] - pop_observed.deaths[i-1];

    double dC = pop_hist[i].getNumReported() - pop_hist[i-1].getNumReported();
    double dR = pop_hist[i].getNumRecoveredReported() - pop_hist[i-1].getNumRecoveredReported();
    double dF = pop_hist[i].getNumFatalReported() - pop_hist[i-1].getNumFatalReported();

    err_sq_total += dC_obs*std::log(dC + eps) - dC;
    err_sq_recov += dR_obs*std::log(dR + eps) - dR;
    err_sq_fatal += dF_obs*std::log(dF + eps) - dF;
  }
  const double denom = (weight_conf + weight_recov + weight_fatal) * nt;
  sub_costs[0] = -scaling * weight_conf * err_sq_total / denom;
  sub_costs[1] = -scaling * weight_recov* err_sq_recov / denom;
  sub_costs[2] = -scaling * weight_fatal* err_sq_fatal / denom;
#endif

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2];

  f_eval_count++;
  return {cost, sub_costs};
}

double
OptimizerPiecewise::getCostGradient(std::vector<double>& grad)
{
  ModelParams params_orig = params; //create a copy of the current parameters

  double norm_sq = 0.0;
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
    norm_sq += grad[j]*grad[j];
  }

  params = params_orig; //Restore current parameters
  copyParam2Vector(params, param_vec);

  return std::sqrt(norm_sq);
}

void
OptimizerPiecewise::updateOptimalSolution(
    const double& cost, const std::array<double,3>& sub_costs, const Vector& param_vec_cur)
{
  std::size_t min_ind = 0;
  for (min_ind = 0; min_ind < cost_min.size(); ++min_ind)
    if (cost < cost_min[min_ind])
      break;

  cost_min.insert(cost_min.begin() + min_ind, cost);
  sub_costs_min.insert(sub_costs_min.begin() + min_ind, sub_costs);
  optimal_param_vec.insert(optimal_param_vec.begin() + min_ind, param_vec_cur);
}

std::vector<ParamBound>
OptimizerPiecewise::getParameterBounds(int nt, int interval_size, bool linear_basis)
{
  int num_nodes = (int)(nt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;

  for (int i = 0; i < num_nodes; ++i)
    bounds[i] = ParamBound(0, 2, delta); //betaN

  for (int j = 1; j <= num_c_params; ++j)
    for (int i = 0; i < num_nodes; ++i)
      bounds[j*num_nodes + i] = ParamBound(0, 1, delta); //(c0 = ce), c1, (c2 = c3)

#if 0
  const int off = (1 + num_c_params)*num_nodes;
  bounds[off  ] = ParamBound(3.0, 3.0, delta); //T_incub0
  bounds[off+1] = ParamBound(2.0, 2.0, delta); //T_incub1
  bounds[off+2] = ParamBound(6.0, 6.0, delta); //T_asympt
  bounds[off+3] = ParamBound(6.0, 6.0, delta); //T_mild
  bounds[off+4] = ParamBound(4.0, 4.0, delta); //T_severe
  bounds[off+5] = ParamBound(10.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.3, 0.3, delta); //f
  bounds[off+7] = ParamBound(0.8, 0.8, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.75, 0.75, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0.02, 0.02, 0.1*delta); //IFR
#endif
  return bounds;
}

void
OptimizerPiecewise::copyParam2Vector(const ModelParams& params, Vector& v)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  constexpr int j0 = OPTIMIZE_C0;
  constexpr int j1 = j0 + OPTIMIZE_C1;
  constexpr int j2 = j1 + OPTIMIZE_C2;

  for (int i = 0; i < num_nodes-1; ++i)
  {
    int ind = i*interval_size;
    v[i] = params.betaN[ind];

    if (OPTIMIZE_C0)
      v[j0*num_nodes + i] = params.c0[ind]; //c0 = ce

    if (OPTIMIZE_C1)
      v[j1*num_nodes + i] = params.c1[ind];

    if (OPTIMIZE_C2)
      v[j2*num_nodes + i] = params.c2[ind]; //c2 = c3
  }
  v[  num_nodes-1] = params.betaN[nt_opt-1];

  if (OPTIMIZE_C0)
    v[j0*num_nodes-1] = params.c0[nt_opt-1];
  if (OPTIMIZE_C1)
    v[j1*num_nodes-1] = params.c1[nt_opt-1];
  if (OPTIMIZE_C2)
    v[j2*num_nodes-1] = params.c2[nt_opt-1];

#if 0
  const int off = (1 + num_c_params)*num_nodes;
  v[off  ] = params.T_incub0;
  v[off+1] = params.T_incub1;
  v[off+2] = params.T_asympt;
  v[off+3] = params.T_mild;
  v[off+4] = params.T_severe;
  v[off+5] = params.T_icu;
  v[off+6] = params.f;
  v[off+7] = params.frac_recover_I1;
  v[off+8] = params.frac_recover_I2;
  v[off+9] = params.IFR;
#endif
}

void
OptimizerPiecewise::copyVector2Param(const Vector& v, ModelParams& params)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  constexpr int j0 = OPTIMIZE_C0;
  constexpr int j1 = j0 + OPTIMIZE_C1;
  constexpr int j2 = j1 + OPTIMIZE_C2;

  for (int i = 0; i < num_nodes-1; ++i)
  {
    const int ind0 = i*interval_size;
    const int ind1 = (i < num_nodes-2) ? (i+1)*interval_size : (nt_opt-1);
    const double L = ind1 - ind0;

    for (int j = 0; j <= L; ++j)
    {
      const double s = linear_basis*(j / L);
      const int ind = ind0 + j;
      params.betaN[ind] = (1-s)*v[i] + s*v[i+1];

      if (OPTIMIZE_C0)
      {
        params.c0[ind]    = (1-s)*v[j0*num_nodes + i] + s*v[j0*num_nodes + i+1];
        params.ce[ind]    = params.c0[ind];
      }
      if (OPTIMIZE_C1)
        params.c1[ind]    = (1-s)*v[j1*num_nodes + i] + s*v[j1*num_nodes + i+1];

      if (OPTIMIZE_C2)
      {
        params.c2[ind]    = (1-s)*v[j2*num_nodes + i] + s*v[j2*num_nodes + i+1];
        params.c3[ind]    = params.c2[ind];
      }
    }
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = params.nt_hist; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = params.betaN[params.nt_hist-1];
    params.ce[i]    = params.ce[params.nt_hist-1];
    params.c0[i]    = params.c0[params.nt_hist-1];
    params.c1[i]    = params.c1[params.nt_hist-1];
    params.c2[i]    = params.c2[params.nt_hist-1];
    params.c3[i]    = params.c3[params.nt_hist-1];
  }
#if 0
  const int off = (1 + num_c_params)*num_nodes;
  params.T_incub0        = v[off  ];
  params.T_incub1        = v[off+1];
  params.T_asympt        = v[off+2];
  params.T_mild          = v[off+3];
  params.T_severe        = v[off+4];
  params.T_icu           = v[off+5];
  params.f               = v[off+6];
  params.frac_recover_I1 = v[off+7];
  params.frac_recover_I2 = v[off+8];
  params.IFR             = v[off+9];
#endif
}

