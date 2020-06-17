#include <iostream>
#include <iomanip>

#include "OptimizerPiecewiseLinear.h"
#include "Simulator.h"
#if 0
OptimizerPiecewiseLinear::OptimizerPiecewiseLinear(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                                 const Vector& quarantine_input, int interval_size_,
                                 double wconf, double wrecov, double wfatal,
                                 int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), interval_size(interval_size_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0, quarantine_input),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  param_bounds = getParameterBounds(nt_opt, interval_size);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}


void
OptimizerPiecewiseLinear::optimizeParametersNLOPT()
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

  cost_min.fill(std::numeric_limits<double>::max());
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

    if (cost_opt < cost_min.back())
    {
      updateOptimalSolution(cost_opt, subcosts_dummy, x);
    }

    randomizeParameters();
    x = param_vec.getDataVector();
  }


  for (int i = max_passes; i < NUM_RESULTS; ++i)
  {
    optimal_param_vec[i] = optimal_param_vec[max_passes-1];
    cost_min[i] = cost_min[max_passes-1];
  }

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerPiecewiseLinear::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerPiecewiseLinear* opt = reinterpret_cast<OptimizerPiecewiseLinear*>(data);

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
OptimizerPiecewiseLinear::getCost()
{
  auto pop_hist = predictModel(params, pop_init);

  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,3> sub_costs = {0.0, 0.0, 0.0};

  double err_sq_total = 0.0, err_sq_recov = 0.0, err_sq_fatal = 0.0;

#if 1
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
#elif 0
  for (int i = 1; i < nt; ++i)
  {
    double trueC = std::max(pop_observed.confirmed[i], 1);
    double trueR = std::max(pop_observed.recovered[i], 1);
    double trueF = std::max(pop_observed.deaths[i], 1);

    double err_total = std::abs(pop_hist[i].getNumReported() - trueC) / trueC;
    double err_recov = std::abs(pop_hist[i].getNumRecoveredReported() - trueR) / trueR;
    double err_fatal = std::abs(pop_hist[i].getNumFatalReported() - trueF) / trueF;

    err_sq_total += err_total*err_total;
    err_sq_recov += err_recov*err_recov;
    err_sq_fatal += err_fatal*err_fatal;
  }
  sub_costs[0] = weight_conf*(err_sq_total)/(3.0*nt);
  sub_costs[1] = weight_recov*(err_sq_recov)/(3.0*nt);
  sub_costs[2] = weight_fatal*(err_sq_fatal)/(3.0*nt);
#else
  const double eps = 1e-10;
  double sq_total = 0.0, sq_recov = 0.0, sq_fatal = 0.0;
  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - pop_observed.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - pop_observed.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - pop_observed.deaths[i]);

    err_sq_total += std::log(std::abs(err_total) + eps);
    err_sq_recov += std::log(std::abs(err_recov) + eps);
    err_sq_fatal += std::log(std::abs(err_fatal) + eps);
  }
  sub_costs[0] = weight_conf*(err_sq_total)/nt;
  sub_costs[1] = weight_recov*(err_sq_recov)/nt;
  sub_costs[2] = weight_fatal*(err_sq_fatal)/nt;
#endif

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2];

  f_eval_count++;
  return {cost, sub_costs};
}

double
OptimizerPiecewiseLinear::getCostGradient(std::vector<double>& grad)
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
OptimizerPiecewiseLinear::updateOptimalSolution(
    const double& cost, const std::array<double,3>& sub_costs, const Vector& param_vec_cur)
{
  int min_ind = 0;
  for (min_ind = 0; min_ind < NUM_RESULTS; ++min_ind)
    if (cost < cost_min[min_ind])
      break;

  for (int i = NUM_RESULTS-1; i > min_ind; i--)
  {
    cost_min[i] = cost_min[i-1];
    sub_costs_min[i] = sub_costs_min[i-1];
    optimal_param_vec[i] = optimal_param_vec[i-1];
  }
  cost_min[min_ind] = cost;
  sub_costs_min[min_ind] = sub_costs;
  optimal_param_vec[min_ind] = param_vec_cur;
}

#if 0
std::vector<ParamBound>
OptimizerPiecewiseLinear::getParameterBounds(int nbasis)
{
  const int m = 5*nbasis + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;
  const double a = 0.1;

  //Bounds for "average" basis coefficients
  bounds[       0] = ParamBound(0, 2, delta); //betaN
  bounds[  nbasis] = ParamBound(0, 1, delta); //c0 = ce
  bounds[2*nbasis] = ParamBound(0, 1, delta); //c1
  bounds[3*nbasis] = ParamBound(0, 1, delta); //c2
  bounds[4*nbasis] = ParamBound(0, 1, delta); //c3

  for (int i = 1; i < nbasis; ++i)
  {
    bounds[           i] = ParamBound(-a, a, delta); //betaN
    bounds[  nbasis + i] = ParamBound(-a, a, delta); //c0 = ce
    bounds[2*nbasis + i] = ParamBound(-a, a, delta); //c1
    bounds[3*nbasis + i] = ParamBound(-a, a, delta); //c2
    bounds[4*nbasis + i] = ParamBound(-a, a, delta); //c3
  }
  const int off = 5*nbasis;
//  bounds[off  ] = ParamBound(3.0, 3.0, delta); //T_incub0
//  bounds[off+1] = ParamBound(2.0, 2.0, delta); //T_incub1
//  bounds[off+2] = ParamBound(6.0, 6.0, delta); //T_asympt
//  bounds[off+3] = ParamBound(6.0, 6.0, delta); //T_mild
//  bounds[off+4] = ParamBound(4.0, 4.0, delta); //T_severe
//  bounds[off+5] = ParamBound(10.0, 10.0, delta); //T_icu
//  bounds[off+6] = ParamBound(0.3, 0.3, delta); //f
//  bounds[off+7] = ParamBound(0.8, 0.8, delta); //frac_recover_I1
//  bounds[off+8] = ParamBound(0.75, 0.75, delta); //frac_recover_I2
//  bounds[off+9] = ParamBound(0.02, 0.02, 0.1*delta); //CFR
//  bounds[off+10] = ParamBound(14.0, 14.0, delta); //T_discharge

//  bounds[off  ] = ParamBound(2.9, 3.1, delta); //T_incub0
//  bounds[off+1] = ParamBound(1.9, 2.1, delta); //T_incub1
//  bounds[off+2] = ParamBound(5.9, 6.1, delta); //T_asympt
//  bounds[off+3] = ParamBound(5.9, 6.1, delta); //T_mild
//  bounds[off+4] = ParamBound(3.9, 4.1, delta); //T_severe
//  bounds[off+5] = ParamBound(9.9, 10.1, delta); //T_icu
//  bounds[off+6] = ParamBound(0.29, 0.31, delta); //f
//  bounds[off+7] = ParamBound(0.5, 0.9, delta); //frac_recover_I1
//  bounds[off+8] = ParamBound(0.5, 0.9, delta); //frac_recover_I2
//  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR
//  bounds[off+10] = ParamBound(13.9, 14.1, delta); //T_discharge

  bounds[off  ] = ParamBound(1.0, 10.0, delta); //T_incub0
  bounds[off+1] = ParamBound(1.0, 10.0, delta); //T_incub1
  bounds[off+2] = ParamBound(1.0, 10.0, delta); //T_asympt
  bounds[off+3] = ParamBound(1.0, 10.0, delta); //T_mild
  bounds[off+4] = ParamBound(1.0, 10.0, delta); //T_severe
  bounds[off+5] = ParamBound(1.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.1, 0.9, delta); //f
  bounds[off+7] = ParamBound(0.1, 0.9, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.1, 0.9, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR
  bounds[off+10] = ParamBound(1.0, 28.0, delta); //T_discharge

  return bounds;
}

void
OptimizerPiecewiseLinear::copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int m = 5*num_basis + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  double avg_beta = 0.0;
  std::array<double, 4> avg_c = {0.0, 0.0, 0.0, 0.0};

  for (int i = 0; i < params.nt_hist; ++i)
  {
    avg_beta += params.betaN[i];
    avg_c[0] += params.c0[i];
    avg_c[1] += params.c1[i];
    avg_c[2] += params.c2[i];
    avg_c[3] += params.c3[i];
  }
  avg_beta /= (double) params.nt_hist;
  v[0] = avg_beta;

  for (std::size_t i = 0; i < avg_c.size(); ++i)
  {
    avg_c[i] /= (double) params.nt_hist;
    v[(i+1)*num_basis] = avg_c[i];
  }

  for (int i = 0; i < 5; ++i)
    for (int j = 1; j < num_basis; ++j)
      v[i*num_basis + j] = 0.0;

  const int off = 5*num_basis;
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
  v[off+10] = params.T_discharge;
}

void
OptimizerPiecewiseLinear::copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 5*num_basis + 11;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  evaluateLegendrePolynomial(num_basis, nt, &v[0], params.betaN);
  evaluateLegendrePolynomial(num_basis, nt, &v[num_basis], params.c0);
  params.ce = params.c0;
  evaluateLegendrePolynomial(num_basis, nt, &v[2*num_basis], params.c1);
  evaluateLegendrePolynomial(num_basis, nt, &v[3*num_basis], params.c2);
  evaluateLegendrePolynomial(num_basis, nt, &v[4*num_basis], params.c3);

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = params.betaN[nt-1];
    params.ce[i]    = params.ce[nt-1];
    params.c0[i]    = params.c0[nt-1];
    params.c1[i]    = params.c1[nt-1];
    params.c2[i]    = params.c2[nt-1];
    params.c3[i]    = params.c3[nt-1];
  }

  const int off = 5*num_basis;
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
  params.T_discharge     = v[off+10];
}
#elif 1
std::vector<ParamBound>
OptimizerPiecewiseLinear::getParameterBounds(int nt, int interval_size)
{
  int num_nodes = (int)(nt/interval_size) + 1;
  const int m = 4*num_nodes + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;
  for (int i = 0; i < num_nodes; ++i)
  {
    bounds[              i] = ParamBound(0, 1, delta); //betaN
    bounds[  num_nodes + i] = ParamBound(0, 1, delta); //c0 = ce
    bounds[2*num_nodes + i] = ParamBound(0, 1, delta); //c1
    bounds[3*num_nodes + i] = ParamBound(1, 1, delta); //c2 = c3
  }
  const int off = 4*num_nodes;
#if 1
  bounds[off  ] = ParamBound(3.0, 3.0, delta); //T_incub0
  bounds[off+1] = ParamBound(2.0, 2.0, delta); //T_incub1
  bounds[off+2] = ParamBound(6.0, 6.0, delta); //T_asympt
  bounds[off+3] = ParamBound(6.0, 6.0, delta); //T_mild
  bounds[off+4] = ParamBound(4.0, 4.0, delta); //T_severe
  bounds[off+5] = ParamBound(10.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.3, 0.3, delta); //f
  bounds[off+7] = ParamBound(0.8, 0.8, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.75, 0.75, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0.02, 0.02, 0.1*delta); //CFR
  bounds[off+10] = ParamBound(14.0, 14.0, delta); //T_discharge
#else
  bounds[off  ] = ParamBound(1.0, 10.0, delta); //T_incub0
  bounds[off+1] = ParamBound(1.0, 10.0, delta); //T_incub1
  bounds[off+2] = ParamBound(1.0, 10.0, delta); //T_asympt
  bounds[off+3] = ParamBound(1.0, 10.0, delta); //T_mild
  bounds[off+4] = ParamBound(1.0, 10.0, delta); //T_severe
  bounds[off+5] = ParamBound(1.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.1, 0.9, delta); //f
  bounds[off+7] = ParamBound(0.1, 0.9, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.1, 0.9, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR
  bounds[off+10] = ParamBound(1.0, 28.0, delta); //T_discharge
#endif

  return bounds;
}

void
OptimizerPiecewiseLinear::copyParam2Vector(const ModelParams& params, Vector& v)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1;
  const int m = 4*num_nodes + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
  {
    int ind = i*interval_size;
    v[              i] = params.betaN[ind];
    v[  num_nodes + i] = params.c0[ind]; //c0 = ce
    v[2*num_nodes + i] = params.c1[ind];
    v[3*num_nodes + i] = params.c2[ind]; //c2 = c3
  }
  v[  num_nodes-1] = params.betaN[nt_opt-1];
  v[2*num_nodes-1] = params.c0[nt_opt-1];
  v[3*num_nodes-1] = params.c1[nt_opt-1];
  v[4*num_nodes-1] = params.c2[nt_opt-1];

  const int off = 4*num_nodes;
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
  v[off+10] = params.T_discharge;
}

void
OptimizerPiecewiseLinear::copyVector2Param(const Vector& v, ModelParams& params)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1;
  const int m = 4*num_nodes + 11;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
  {
    const int ind0 = i*interval_size;
    const int ind1 = (i < num_nodes-2) ? (i+1)*interval_size : (nt_opt-1);
    const double L = ind1 - ind0;

    for (int j = 0; j <= L; ++j)
    {
      const double s = j / L;
      const int ind = ind0 + j;
      params.betaN[ind] = (1-s)*v[              i] + s*v[              i+1];
      params.c0[ind]    = (1-s)*v[  num_nodes + i] + s*v[  num_nodes + i+1];
      params.ce[ind]    = params.c0[ind];
      params.c1[ind]    = (1-s)*v[2*num_nodes + i] + s*v[2*num_nodes + i+1];
      params.c2[ind]    = (1-s)*v[3*num_nodes + i] + s*v[3*num_nodes + i+1];
      params.c3[ind]    = params.c2[ind];
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

  const int off = 4*num_nodes;
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
  params.T_discharge     = v[off+10];
}

#else
std::vector<ParamBound>
OptimizerPiecewiseLinear::getParameterBounds(int nt, int interval_size)
{
  int num_nodes = (int)(nt/interval_size) + 1;
  const int m = num_nodes + 3 + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;
  for (int i = 0; i < num_nodes; ++i)
    bounds[i] = ParamBound(0, 1, delta); //betaN
  bounds[num_nodes  ] = ParamBound(0, 1, delta); //c0 = ce
  bounds[num_nodes+1] = ParamBound(0, 1, delta); //c1
  bounds[num_nodes+2] = ParamBound(1, 1, delta); //c2 = c3

  const int off = num_nodes + 3;
#if 1
  bounds[off  ] = ParamBound(3.0, 3.0, delta); //T_incub0
  bounds[off+1] = ParamBound(2.0, 2.0, delta); //T_incub1
  bounds[off+2] = ParamBound(6.0, 6.0, delta); //T_asympt
  bounds[off+3] = ParamBound(6.0, 6.0, delta); //T_mild
  bounds[off+4] = ParamBound(4.0, 4.0, delta); //T_severe
  bounds[off+5] = ParamBound(10.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.3, 0.3, delta); //f
  bounds[off+7] = ParamBound(0.8, 0.8, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.75, 0.75, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0.02, 0.02, 0.1*delta); //CFR
  bounds[off+10] = ParamBound(14.0, 14.0, delta); //T_discharge
#else
  bounds[off  ] = ParamBound(1.0, 10.0, delta); //T_incub0
  bounds[off+1] = ParamBound(1.0, 10.0, delta); //T_incub1
  bounds[off+2] = ParamBound(1.0, 10.0, delta); //T_asympt
  bounds[off+3] = ParamBound(1.0, 10.0, delta); //T_mild
  bounds[off+4] = ParamBound(1.0, 10.0, delta); //T_severe
  bounds[off+5] = ParamBound(1.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.1, 0.9, delta); //f
  bounds[off+7] = ParamBound(0.1, 0.9, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.1, 0.9, delta); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR
  bounds[off+10] = ParamBound(1.0, 28.0, delta); //T_discharge
#endif

  return bounds;
}

void
OptimizerPiecewiseLinear::copyParam2Vector(const ModelParams& params, Vector& v)
{
  int num_nodes = (int)(nt_opt/num_basis) + 1;
  const int m = num_nodes + 3 + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
    v[i] = params.betaN[i*num_basis];
  v[num_nodes-1] = params.betaN[nt_opt-1];

  v[num_nodes  ] = params.c0[0]; //c0 = ce
  v[num_nodes+1] = params.c1[0]; //c1
  v[num_nodes+2] = params.c2[0]; //c2 = c3

  const int off = num_nodes + 3;
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
  v[off+10] = params.T_discharge;
}

void
OptimizerPiecewiseLinear::copyVector2Param(const Vector& v, ModelParams& params)
{
  int num_nodes = (int)(nt_opt/num_basis) + 1;
  const int m = num_nodes + 3 + 11;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
  {
    const int ind0 = i*num_basis;
    const int ind1 = (i < num_nodes-2) ? (i+1)*num_basis : (nt_opt-1);
    const double L = ind1 - ind0;

    for (int j = 0; j <= L; ++j)
    {
      const double s = j / L;
      params.betaN[ind0 + j] = (1-s)*v[i] + s*v[i+1];
    }
  }

  for (int i = 0; i < params.nt_hist; ++i)
  {
    params.ce[i] = v[num_nodes];
    params.c0[i] = params.ce[i];
    params.c1[i] = v[num_nodes+1];
    params.c2[i] = v[num_nodes+2];
    params.c3[i] = params.c2[i];
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

  const int off = num_nodes + 3;
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
  params.T_discharge     = v[off+10];
}
#endif

void
OptimizerPiecewiseLinear::evaluateLegendrePolynomial(const int nbasis, const int nt, const double* coeff, Vector& params)
{
  if (nbasis == 1)
  {
    for (int i = 0; i < nt; ++i)
      params[i] = coeff[0];
  }
  else if (nbasis == 2)
  {
    for (int i = 0; i < nt; ++i)
    {
      double x = 2.0*(i / (double)(nt - 1)) - 1.0; //in range [-1, 1]
      double phi0 = 1.0;
      double phi1 = x;
      params[i] = coeff[0]*phi0 + coeff[1]*phi1;
    }
  }
  else if (nbasis == 3)
  {
    for (int i = 0; i < nt; ++i)
    {
      double x = 2.0*(i / (double)(nt - 1)) - 1.0; //in range [-1, 1]
      double phi0 = 1.0;
      double phi1 = x;
      double phi2 = 1.5*x*x - 0.5;
      params[i] = coeff[0]*phi0 + coeff[1]*phi1 + coeff[2]*phi2;
    }
  }
  else
    throwError("evaluateLegendrePolynomial - Unsupported basis number!");
}
#endif
