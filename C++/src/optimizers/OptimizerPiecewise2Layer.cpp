#include <iostream>
#include <iomanip>

#include "OptimizerPiecewise2Layer.h"
#include "models/Simulator.h"
#include "linearalgebra/SurrealS.h"

constexpr std::array<int,9> OptimizerPiecewise2Layer::IFR_SEG_STARTS;

OptimizerPiecewise2Layer::OptimizerPiecewise2Layer(
    const ObservedPopulation& pop_observed_, const Population2Layer<double>& pop_init_,
    const Vector<double>& testing_data, 
    const Vector<double>& vaccination_data,
    int interval_size_, bool linear_basis_,
    double wconf, double wrecov, double wfatal,
    int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), interval_size(interval_size_), linear_basis(linear_basis_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0, testing_data, vaccination_data),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  param_bounds = getParameterBounds(nt_opt, interval_size, linear_basis);
//  param_vec.resize(param_bounds.size());
//  copyParam2Vector(params, param_vec);
}

void
OptimizerPiecewise2Layer::optimizeParametersNLOPT()
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
OptimizerPiecewise2Layer::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerPiecewise2Layer* opt = reinterpret_cast<OptimizerPiecewise2Layer*>(data);
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
OptimizerPiecewise2Layer::getCost(const ModelParams2Layer<T>& params_tmp)
{
  using namespace std;

  auto pop_hist = predictModel(params_tmp, pop_init);

  const int nt = params_tmp.nt_hist + params_tmp.nt_pred;

  std::array<T,3> sub_costs = {0.0, 0.0, 0.0};
  T err_sq_total = 0.0, err_sq_fatal = 0.0;

#if 1 //L2
  double sq_total = 0.0, sq_fatal = 0.0;
  double obs_confirmed, obs_deaths;
  T err_total, err_fatal;
  for (int i = 0; i < nt; ++i)
  {
    obs_confirmed = pop_observed.confirmed[i];
    obs_deaths    = pop_observed.deaths[i];

    err_total = (pop_hist[i].getNumReported() - obs_confirmed);
    err_fatal = (pop_hist[i].getNumFatalReported() - obs_deaths);

    err_sq_total += err_total*err_total;
    err_sq_fatal += err_fatal*err_fatal;

    sq_total += obs_confirmed*obs_confirmed;
    sq_fatal += obs_deaths   *obs_deaths;
  }
  sub_costs[0] = weight_conf*(err_sq_total/sq_total);
  //sub_costs[1] = weight_recov*(err_sq_recov/sq_recov);
  sub_costs[2] = weight_fatal*(err_sq_fatal/sq_fatal);
#elif 0 //L2 of log
  for (int i = 1; i < nt; ++i)
  {
    T logC_obs = log(max(pop_observed.confirmed[i], 1));
    T logF_obs = log(max(pop_observed.deaths[i], 1));

    T logC = log(max(pop_hist[i].getNumReported(), 1.0));
    T logF = log(max(pop_hist[i].getNumFatalReported(), 1.0));

    err_sq_total += (logC - logC_obs)*(logC - logC_obs);
    err_sq_fatal += (logF - logF_obs)*(logF - logF_obs);
  }
  const double denom = (weight_conf + weight_recov + weight_fatal) * nt;
  sub_costs[0] = weight_conf * err_sq_total / denom;
  // sub_costs[1] = weight_recov* err_sq_recov / denom;
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
OptimizerPiecewise2Layer::getCostGradientFD(
  const ModelParams2Layer<double>& params_orig, std::vector<double>& grad)
{
  ModelParams2Layer<double> params_tmp = params_orig; //create a copy of the current parameters

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
OptimizerPiecewise2Layer::getCostGradientAD(
  const ModelParams2Layer<double>& params_orig, std::vector<double>& grad)
{
  ModelParams2Layer<SurrealS<1,double>> params_tmp(params_orig); //create a copy of the current parameters
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
OptimizerPiecewise2Layer::updateOptimalSolution(
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
OptimizerPiecewise2Layer::getParameterBounds(int nt, int interval_size, bool linear_basis)
{
  int num_nodes = (int)(nt/interval_size) + 1*linear_basis;
  constexpr int num_dynamic_params = 1;
  //[beta values, IFR segment values, T_incub, T_recov, beta_test_scaling, beta_vac_scaling, vaccine_alpha, IFR_vac_scaling]
  const int m = num_dynamic_params*num_nodes + IFR_SEG_STARTS.size() + 6;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;

  for (int i = 0; i < num_nodes; ++i)
    bounds[i] = ParamBound(0, 1.0, delta); //beta

  int off = num_dynamic_params*num_nodes;
  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
    bounds[off + i] = ParamBound(0, 0.02, delta); //IFR in [0%, 2%]

  off += IFR_SEG_STARTS.size();
  bounds[off++] = ParamBound(2.0, 10.0, delta); //T_incub
  bounds[off++] = ParamBound(2.0, 14.0, delta); //T_recov
  bounds[off++] = ParamBound(0.0, 1.0, delta); //beta_test_scaling
  bounds[off++] = ParamBound(0.0, 1.0, delta); //beta_vac_scaling
  bounds[off++] = ParamBound(0.0, 1.0, delta); //vaccine_alpha
  bounds[off++] = ParamBound(0.0, 1.0, delta); //IFR_vac_scaling

  return bounds;
}

template<class T>
void
OptimizerPiecewise2Layer::copyParam2Vector(const ModelParams2Layer<T>& params, std::vector<T>& v)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_dynamic_params = 1;
  //[beta values, IFR segment values, T_incub, T_recov, beta_test_scaling, beta_vac_scaling, vaccine_alpha, IFR_vac_scaling]
  const int m = num_dynamic_params*num_nodes + IFR_SEG_STARTS.size() + 6;
  if ((int) v.size() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
  {
    int ind = i*interval_size;
    v[i] = params.beta[ind];
  }
  v[num_nodes-1] = params.beta[nt_opt-1];

  int off = num_dynamic_params*num_nodes;
  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
    v[off + i] = params.IFR[IFR_SEG_STARTS[i]];

  off += IFR_SEG_STARTS.size();
  v[off++] = params.T_incub;
  v[off++] = params.T_recov;
  v[off++] = params.beta_test_scaling;
  v[off++] = params.beta_vac_scaling;
  v[off++] = params.vaccine_alpha;
  v[off++] = params.IFR_vac_scaling;
}

template<class T>
void
OptimizerPiecewise2Layer::copyVector2Param(const std::vector<T>& v, ModelParams2Layer<T>& params)
{
  int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
  constexpr int num_dynamic_params = 1;
  //[beta values, IFR segment values, T_incub, T_recov, beta_test_scaling, beta_vac_scaling, vaccine_alpha, IFR_vac_scaling]
  const int m = num_dynamic_params*num_nodes + IFR_SEG_STARTS.size() + 6;
  if ((int)v.size() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

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
        params.beta[ind] = (1-s)*v[i] + s*v[i+1];
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
        params.beta[ind] = v[i];
    }
  }

  int off = num_dynamic_params*num_nodes;
  for (std::size_t i = 0; i < IFR_SEG_STARTS.size(); ++i)
  {
    const int ind1 = ((int)i < (int)IFR_SEG_STARTS.size()-1) ? IFR_SEG_STARTS[i+1] : nt_opt;
    for (int ind = IFR_SEG_STARTS[i]; ind < ind1; ++ind)
      params.IFR[ind] = v[off + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = params.nt_hist; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.beta[i] = params.beta[params.nt_hist-1];
    params.IFR[i]  = params.IFR[params.nt_hist-1];
  }

  off += IFR_SEG_STARTS.size();
  params.T_incub           = v[off++];
  params.T_recov           = v[off++];
  params.beta_test_scaling = v[off++];
  params.beta_vac_scaling  = v[off++];
  params.vaccine_alpha     = v[off++];
  params.IFR_vac_scaling   = v[off++];
}
