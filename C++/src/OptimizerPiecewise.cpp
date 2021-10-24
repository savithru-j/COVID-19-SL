#include <iostream>
#include <iomanip>

#include "OptimizerPiecewise.h"
#include "Simulator.h"
#include "SurrealS.h"

constexpr std::array<int,3> OptimizerPiecewise::IFR_SEG_STARTS;

OptimizerPiecewise::OptimizerPiecewise(
    const ObservedPopulation& pop_observed_, const Population<double>& pop_init_,
    const Vector<double>& quarantine_input, const Vector<double>& vaccination_data,
    int interval_size_, bool linear_basis_,
    double wconf, double wrecov, double wfatal,
    int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), interval_size(interval_size_), linear_basis(linear_basis_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0, quarantine_input, vaccination_data),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  param_bounds = getParameterBounds(nt_opt, interval_size, linear_basis);
//  param_vec.resize(param_bounds.size());
//  copyParam2Vector(params, param_vec);
}

void
OptimizerPiecewise::optimizeParametersNLOPT()
{
  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;
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

  std::vector<double> x(nDim());
  copyParam2Vector(params, x);

  std::cout << std::scientific << std::setprecision(5);

  for (int pass = 0; pass < max_passes; ++pass)
  {
    std::cout << std::endl << "Pass " << (pass+1) << "/" << max_passes << " ----------------" << std::endl;

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

    randomizeParameters(x);
  }

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerPiecewise::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerPiecewise* opt = reinterpret_cast<OptimizerPiecewise*>(data);
  opt->copyVector2Param(x, opt->params);

  auto cost = opt->getCost(opt->params);

  double grad_norm = 0.0;
  if (!grad.empty())
    grad_norm = opt->getCostGradientAD(opt->params, grad);

  std::cout << "  Iter: " << opt->nlopt_iter << ", Cost: " << cost.first
            << " = [" << cost.second << "], Grad-norm: " << grad_norm << std::endl;
  opt->nlopt_iter++;

  return cost.first;
}

template<class T>
std::pair<T, std::array<T, 3>>
OptimizerPiecewise::getCost(const ModelParams<T>& params_tmp)
{
  using namespace std;

  auto pop_hist = predictModel(params_tmp, pop_init);

  const int nt = params_tmp.nt_hist + params_tmp.nt_pred;

  std::array<T,3> sub_costs = {0.0, 0.0, 0.0};
  double obs_confirmed, obs_recovered, obs_deaths;
  T err_total, err_recov, err_fatal;
  T err_sq_total = 0.0, err_sq_recov = 0.0, err_sq_fatal = 0.0;

#if 1 //L2
  double sq_total = 0.0, sq_recov = 0.0, sq_fatal = 0.0;
  for (int i = 0; i < nt; ++i)
  {
    obs_confirmed = pop_observed.confirmed[i];
    obs_recovered = pop_observed.recovered[i];
    obs_deaths    = pop_observed.deaths[i];

    err_total = (pop_hist[i].getNumReported() - obs_confirmed);
    err_recov = (pop_hist[i].getNumRecoveredReported() - obs_recovered);
    err_fatal = (pop_hist[i].getNumFatalReported() - obs_deaths);

    err_sq_total += err_total*err_total;
    err_sq_recov += err_recov*err_recov;
    err_sq_fatal += err_fatal*err_fatal;

    sq_total += obs_confirmed*obs_confirmed;
    sq_recov += obs_recovered*obs_recovered;
    sq_fatal += obs_deaths   *obs_deaths;
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
#elif 0 //MLE
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

  const T cost = sub_costs[0] + sub_costs[1] + sub_costs[2];

  f_eval_count++;
  return {cost, sub_costs};
}

double
OptimizerPiecewise::getCostGradientFD(const ModelParams<double>& params_orig, std::vector<double>& grad)
{
  ModelParams<double> params_tmp = params_orig; //create a copy of the current parameters

  Vector<double> param_vec(nDim());
  copyParam2Vector(params_tmp, param_vec);

  double norm_sq = 0.0;
  for (int j = 0; j < param_vec.m(); ++j)
  {
    const double delta = param_bounds[j].step;

    //Compute finite difference
    param_vec[j] += delta;
    copyVector2Param(param_vec, params_tmp);
    double fp = getCost(params_tmp).first;

    param_vec[j] -= 2*delta;
    copyVector2Param(param_vec, params_tmp);
    double fm = getCost(params_tmp).first;

    param_vec[j] += delta;

    grad[j] = (fp - fm)/(2*delta);
    norm_sq += grad[j]*grad[j];
  }

  return std::sqrt(norm_sq);
}

double
OptimizerPiecewise::getCostGradientAD(const ModelParams<double>& params_orig, std::vector<double>& grad)
{
  ModelParams<SurrealS<1,double>> params_tmp(params_orig); //create a copy of the current parameters
  Vector<SurrealS<1,double>> param_vec(nDim());
  copyParam2Vector(params_tmp, param_vec);

  double norm_sq = 0.0;
  for (int j = 0; j < param_vec.m(); ++j)
  {
    param_vec[j].deriv(0) = 1.0; //derivative w.r.t parameter j
    copyVector2Param(param_vec, params_tmp);
    auto f = getCost(params_tmp).first;
    grad[j] = f.deriv(0);
    norm_sq += grad[j]*grad[j];
    param_vec[j].deriv(0) = 0.0; //reset derivative
  }
  return std::sqrt(norm_sq);
}

void
OptimizerPiecewise::updateOptimalSolution(
    const double& cost, const std::array<double,3>& sub_costs, const Vector<double>& param_vec_cur)
{
  std::size_t min_ind = 0;
  for (min_ind = 0; min_ind < cost_min.size(); ++min_ind)
    if (cost < cost_min[min_ind])
      break;

  cost_min.insert(cost_min.begin() + min_ind, cost);
  sub_costs_min.insert(sub_costs_min.begin() + min_ind, sub_costs);
  optimal_param_vec.insert(optimal_param_vec.begin() + min_ind, param_vec_cur);
}

Vector<ParamBound>
OptimizerPiecewise::getParameterBounds(int nt, int interval_size, bool linear_basis)
{
  int num_nodes = (int)(nt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes + IFR_SEG_STARTS.size() + 1; //[beta values, c_params values, IFR segment values, vaccine_eff]
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;

  for (int i = 0; i < num_nodes; ++i)
    bounds[i] = ParamBound(0, 2, delta); //betaN

  for (int j = 1; j <= num_c_params; ++j)
    for (int i = 0; i < num_nodes; ++i)
      bounds[j*num_nodes + i] = ParamBound(0, 1, delta); //(c0 = ce), c1, (c2 = c3)

  int off = (1 + num_c_params)*num_nodes;

  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
    bounds[off + i] = ParamBound(0, 0.02, delta); //IFR in [0%, 2%]

  off += IFR_SEG_STARTS.size();
  bounds[off] = ParamBound(0.0, 1.0, delta); //vaccine_eff

#if 0
  off += 1;
  bounds[off  ] = ParamBound(3.0, 3.0, delta); //T_incub0
  bounds[off+1] = ParamBound(2.0, 2.0, delta); //T_incub1
  bounds[off+2] = ParamBound(6.0, 6.0, delta); //T_asympt
  bounds[off+3] = ParamBound(6.0, 6.0, delta); //T_mild
  bounds[off+4] = ParamBound(4.0, 4.0, delta); //T_severe
  bounds[off+5] = ParamBound(10.0, 10.0, delta); //T_icu
  bounds[off+6] = ParamBound(0.3, 0.3, delta); //f
  bounds[off+7] = ParamBound(0.8, 0.8, delta); //frac_recover_I1
  bounds[off+8] = ParamBound(0.75, 0.75, delta); //frac_recover_I2
#endif
  return bounds;
}

template<class T>
void
OptimizerPiecewise::copyParam2Vector(const ModelParams<T>& params, std::vector<T>& v)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes + IFR_SEG_STARTS.size() + 1; //[beta values, c_params values, IFR segment values, vaccine_eff]
  if ((int) v.size() != m)
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
    v[(j0+1)*num_nodes-1] = params.c0[nt_opt-1];
  if (OPTIMIZE_C1)
    v[(j1+1)*num_nodes-1] = params.c1[nt_opt-1];
  if (OPTIMIZE_C2)
    v[(j2+1)*num_nodes-1] = params.c2[nt_opt-1];

  int off = (1 + num_c_params)*num_nodes;

  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
    v[off + i] = params.IFR[IFR_SEG_STARTS[i]];

  off += IFR_SEG_STARTS.size();
  v[off] = params.vaccine_eff;

#if 0
  off += 1;
  v[off  ] = params.T_incub0;
  v[off+1] = params.T_incub1;
  v[off+2] = params.T_asympt;
  v[off+3] = params.T_mild;
  v[off+4] = params.T_severe;
  v[off+5] = params.T_icu;
  v[off+6] = params.f;
  v[off+7] = params.frac_recover_I1;
  v[off+8] = params.frac_recover_I2;
#endif
}

template<class T>
void
OptimizerPiecewise::copyVector2Param(const std::vector<T>& v, ModelParams<T>& params)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
  const int m = (1 + num_c_params)*num_nodes + IFR_SEG_STARTS.size() + 1; //[beta values, c_params values, IFR segment values, vaccine_eff]
  if ((int)v.size() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  constexpr int j0 = OPTIMIZE_C0;
  constexpr int j1 = j0 + OPTIMIZE_C1;
  constexpr int j2 = j1 + OPTIMIZE_C2;

  if (linear_basis)
  {
    for (int i = 0; i < num_nodes-1; ++i)
    {
      const int ind0 = i*interval_size;
      const int ind1 = (i < num_nodes-2) ? (i+1)*interval_size : (nt_opt-1);
      const double L = ind1 - ind0;

      for (int j = 0; j <= L; ++j)
      {
        const double s = (double) j / L;
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
  }
  else
  {
    for (int i = 0; i < num_nodes; ++i)
    {
      const int ind0 = i*interval_size;
      const int ind1 = (i < num_nodes-1) ? (i+1)*interval_size : nt_opt;
      for (int ind = ind0; ind < ind1; ++ind)
      {
        params.betaN[ind] = v[i];

        if (OPTIMIZE_C0)
        {
          params.c0[ind] = v[j0*num_nodes + i];
          params.ce[ind] = params.c0[ind];
        }
        if (OPTIMIZE_C1)
          params.c1[ind] = v[j1*num_nodes + i];

        if (OPTIMIZE_C2)
        {
          params.c2[ind] = v[j2*num_nodes + i];
          params.c3[ind] = params.c2[ind];
        }
      }
    }
  }

  int off = (1 + num_c_params)*num_nodes;

  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
  {
    const int ind1 = ((int)i < (int)IFR_SEG_STARTS.size()-1) ? IFR_SEG_STARTS[i+1] : nt_opt;
    for (int ind = IFR_SEG_STARTS[i]; ind < ind1; ++ind)
      params.IFR[ind] = v[off + i];
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
    params.IFR[i]   = params.IFR[params.nt_hist-1];
  }

  off += IFR_SEG_STARTS.size();
  params.vaccine_eff = v[off];

#if 0
  off += 1;
  params.T_incub0        = v[off  ];
  params.T_incub1        = v[off+1];
  params.T_asympt        = v[off+2];
  params.T_mild          = v[off+3];
  params.T_severe        = v[off+4];
  params.T_icu           = v[off+5];
  params.f               = v[off+6];
  params.frac_recover_I1 = v[off+7];
  params.frac_recover_I2 = v[off+8];
#endif
}
