#include <iostream>
#include <iomanip>
#include <random>

#include "Optimizer.h"
#include "Simulator.h"

Optimizer::Optimizer(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                     double wconf, double wrecov, double wfatal, double wreg,
                     int max_iter_per_pass_, int max_passes_) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal), weight_reg(wreg),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    nt_opt(pop_observed.getNumDays() - t_buffer),
    params(nt_opt, t_buffer)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  if (weight_reg != 0.0)
  {
//    reg_matrix = getHaarMatrix(nt_opt);
    reg_matrix = getDCTMatrix(nt_opt);
  }

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

  cost_rel_min.fill(std::numeric_limits<double>::max());
  std::array<double,6> subcosts_dummy;

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

    if (cost_opt < cost_rel_min.back())
    {
      updateOptimalSolution(cost_opt, subcosts_dummy, x);
    }

    randomizeParameters(1.0);
    x = param_vec.getDataVector();
  }


  for (int i = max_passes; i < NUM_RESULTS; ++i)
    optimal_param_vec[i] = optimal_param_vec[max_passes-1];

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_rel_min << std::endl;

  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
Optimizer::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  Optimizer* opt = reinterpret_cast<Optimizer*>(data);

  opt->param_vec = x;
  copyVector2Param(opt->param_vec, opt->params);

  auto cost = opt->getCost();

  double grad_norm = 0.0;
  if (!grad.empty())
    grad_norm = opt->getCostGradient(grad);

  std::cout << "  Iter: " << opt->nlopt_iter << ", Cost: " << cost.first
            << " = [" << cost.second << "], Grad-norm: " << grad_norm << std::endl;
  opt->nlopt_iter++;

  return cost.first;
}


std::pair<double, std::array<double, 6>>
Optimizer::getCost()
{
  auto pop_hist = predictModel(params, pop_init);

  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,6> sub_costs = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double err_sq_total = 0.0, err_sq_recov = 0.0, err_sq_fatal = 0.0;
  double sq_total = 0.0, sq_recov = 0.0, sq_fatal = 0.0;

#if 1
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
  sub_costs[0] = weight_conf*(err_sq_total/sq_total) + weight_recov*(err_sq_recov/sq_recov) + weight_fatal*(err_sq_fatal/sq_fatal);
#else
  const double eps = 1e-10;
  for (int i = 1; i < nt; ++i)
  {
    double logC = std::log(pop_observed.confirmed[i] + eps);
    double logR = std::log(pop_observed.recovered[i] + eps);
    double logF = std::log(pop_observed.deaths[i] + eps);

    double err_total = (std::log(pop_hist[i].getNumReported() + eps) - logC);
    double err_recov = (std::log(pop_hist[i].getNumRecoveredReported() + eps) - logR);
    double err_fatal = (std::log(pop_hist[i].getNumFatalReported() + eps) - logF);

    err_sq_total += err_total*err_total;
    err_sq_recov += err_recov*err_recov;
    err_sq_fatal += err_fatal*err_fatal;

    sq_total += logC*logC;
    sq_recov += logR*logR;
    sq_fatal += logF*logF;
  }
  sub_costs[0] = weight_conf*(err_sq_total/sq_total) + weight_recov*(err_sq_recov/sq_recov) + weight_fatal*(err_sq_fatal/sq_fatal);
#endif

  //Add regularization terms
  if (weight_reg != 0.0)
  {
    for (int i = 0; i < reg_matrix.m(); ++i)
    {
      double coeff_betaN = 0.0;
      double coeff_c0 = 0.0;
      double coeff_c1 = 0.0;
      double coeff_c2 = 0.0;
      double coeff_c3 = 0.0;

      for (int j = 0; j < reg_matrix.n(); ++j)
      {
        if (j < params.betaN.m())
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN(j);
          coeff_c0 += reg_matrix(i,j) * params.c0(j);
          coeff_c1 += reg_matrix(i,j) * params.c1(j);
          coeff_c2 += reg_matrix(i,j) * params.c2(j);
          coeff_c3 += reg_matrix(i,j) * params.c3(j);
        }
        else
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN.back();
          coeff_c0 += reg_matrix(i,j) * params.c0.back();
          coeff_c1 += reg_matrix(i,j) * params.c1.back();
          coeff_c2 += reg_matrix(i,j) * params.c2.back();
          coeff_c3 += reg_matrix(i,j) * params.c3.back();
        }
      }

      //Take L1-norms of coefficient vectors
      sub_costs[1] += std::abs(coeff_betaN);
      sub_costs[2] += std::abs(coeff_c0);
      sub_costs[3] += std::abs(coeff_c1);
      sub_costs[4] += std::abs(coeff_c2);
      sub_costs[5] += std::abs(coeff_c3);
    }

    for (int i = 1; i < 6; ++i)
      sub_costs[i] *= weight_reg;
  } //if-regularize

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2] + sub_costs[3] + sub_costs[4] + sub_costs[5];

  f_eval_count++;
  return {cost, sub_costs};
}

double
Optimizer::getCostGradient(std::vector<double>& grad)
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
Optimizer::updateOptimalSolution(
    const double& cost_rel, const std::array<double,6>& sub_costs, const Vector& param_vec_cur)
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

void
Optimizer::copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int nt = params.nt_hist;
  const int m = 5*nt + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");
//  if (v.m() != m)
//    v.resize(m);

  for (int i = 0; i < nt; ++i)
  {
    v[       i] = params.betaN[i];
    v[  nt + i] = params.c0[i];
    v[2*nt + i] = params.c1[i];
    v[3*nt + i] = params.c2[i];
    v[4*nt + i] = params.c3[i];
  }
  const int off = 5*nt;
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
Optimizer::copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 5*nt + 11;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[       i];
    params.ce[i]    = v[  nt + i];
    params.c0[i]    = v[  nt + i]; //ce = c0
    params.c1[i]    = v[2*nt + i];
    params.c2[i]    = v[3*nt + i];
    params.c3[i]    = v[4*nt + i];
  }

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

  const int off = 5*nt;
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

std::vector<ParamBound>
Optimizer::getParameterBounds(int nt)
{
  const int m = 5*nt + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;

  for (int i = 0; i < nt; ++i)
  {
    bounds[       i] = ParamBound(0, 2, delta); //betaN
    bounds[  nt + i] = ParamBound(0, 1, delta); //c0 = ce
    bounds[2*nt + i] = ParamBound(0, 1, delta); //c1
    bounds[3*nt + i] = ParamBound(0, 1, delta); //c2
    bounds[4*nt + i] = ParamBound(0, 1, delta); //c3
  }
  const int off = 5*nt;
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

//  bounds[off  ] = ParamBound(1.0, 10.0, delta); //T_incub0
//  bounds[off+1] = ParamBound(1.0, 10.0, delta); //T_incub1
//  bounds[off+2] = ParamBound(1.0, 10.0, delta); //T_asympt
//  bounds[off+3] = ParamBound(1.0, 10.0, delta); //T_mild
//  bounds[off+4] = ParamBound(1.0, 10.0, delta); //T_severe
//  bounds[off+5] = ParamBound(1.0, 10.0, delta); //T_icu
//  bounds[off+6] = ParamBound(0.1, 0.9, delta); //f
//  bounds[off+7] = ParamBound(0.1, 0.9, delta); //frac_recover_I1
//  bounds[off+8] = ParamBound(0.1, 0.9, delta); //frac_recover_I2
//  bounds[off+9] = ParamBound(0, 0.02, 0.1*delta); //CFR
//  bounds[off+10] = ParamBound(1.0, 28.0, delta); //T_discharge

  return bounds;
}

OptimizerLowDim::OptimizerLowDim(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                                 int num_basis_, double wconf, double wrecov, double wfatal,
                                 int max_iter_per_pass_, int max_passes_) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), num_basis(num_basis_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  param_bounds = getParameterBounds(nt_opt, num_basis);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}


void
OptimizerLowDim::optimizeParametersNLOPT()
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
    optimal_param_vec[i] = optimal_param_vec[max_passes-1];

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerLowDim::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerLowDim* opt = reinterpret_cast<OptimizerLowDim*>(data);

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
OptimizerLowDim::getCost()
{
  auto pop_hist = predictModel(params, pop_init);

  const int nt = params.nt_hist + params.nt_pred;

  std::array<double,3> sub_costs = {0.0, 0.0, 0.0};

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

    sq_total += (double)pop_observed.confirmed[i]*(double)pop_observed.confirmed[i];
    sq_recov += (double)pop_observed.recovered[i]*(double)pop_observed.recovered[i];
    sq_fatal += (double)pop_observed.deaths[i]   *(double)pop_observed.deaths[i];
  }
  sub_costs[0] = weight_conf*(err_sq_total/sq_total);
  sub_costs[1] = weight_recov*(err_sq_recov/sq_recov);
  sub_costs[2] = weight_fatal*(err_sq_fatal/sq_fatal);

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2];

  f_eval_count++;
  return {cost, sub_costs};
}

double
OptimizerLowDim::getCostGradient(std::vector<double>& grad)
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
OptimizerLowDim::updateOptimalSolution(
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
OptimizerLowDim::getParameterBounds(int nbasis)
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
OptimizerLowDim::copyParam2Vector(const ModelParams& params, Vector& v)
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
OptimizerLowDim::copyVector2Param(const Vector& v, ModelParams& params)
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
#endif

std::vector<ParamBound>
OptimizerLowDim::getParameterBounds(int nt, int interval_size)
{
  int num_nodes = (int)(nt/interval_size) + 1;
  const int m = 5*num_nodes + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;
  for (int i = 0; i < num_nodes; ++i)
  {
    bounds[              i] = ParamBound(0, 2, delta); //betaN
    bounds[  num_nodes + i] = ParamBound(0, 1, delta); //c0 = ce
    bounds[2*num_nodes + i] = ParamBound(0, 1, delta); //c1
    bounds[3*num_nodes + i] = ParamBound(0, 1, delta); //c2
    bounds[4*num_nodes + i] = ParamBound(0, 1, delta); //c3
  }
  const int off = 5*num_nodes;
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
OptimizerLowDim::copyParam2Vector(const ModelParams& params, Vector& v)
{
  int num_nodes = (int)(nt_opt/num_basis) + 1;
  const int m = 5*num_nodes + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  for (int i = 0; i < num_nodes-1; ++i)
  {
    int ind = i*num_basis;
    v[              i] = params.betaN[ind];
    v[  num_nodes + i] = params.c0[ind];
    v[2*num_nodes + i] = params.c1[ind];
    v[3*num_nodes + i] = params.c2[ind];
    v[4*num_nodes + i] = params.c3[ind];
  }
  v[  num_nodes-1] = params.betaN[nt_opt-1];
  v[2*num_nodes-1] = params.c0[nt_opt-1];
  v[3*num_nodes-1] = params.c1[nt_opt-1];
  v[4*num_nodes-1] = params.c2[nt_opt-1];
  v[5*num_nodes-1] = params.c3[nt_opt-1];

  const int off = 5*num_nodes;
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
OptimizerLowDim::copyVector2Param(const Vector& v, ModelParams& params)
{
  int num_nodes = (int)(nt_opt/num_basis) + 1;
  const int m = 5*num_nodes + 11;
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
      const int ind = ind0 + j;
      params.betaN[ind] = (1-s)*v[              i] + s*v[              i+1];
      params.c0[ind]    = (1-s)*v[  num_nodes + i] + s*v[  num_nodes + i+1];
      params.c1[ind]    = (1-s)*v[2*num_nodes + i] + s*v[2*num_nodes + i+1];
      params.c2[ind]    = (1-s)*v[3*num_nodes + i] + s*v[3*num_nodes + i+1];
      params.c3[ind]    = (1-s)*v[4*num_nodes + i] + s*v[4*num_nodes + i+1];
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

  const int off = 5*num_nodes;
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


void
OptimizerLowDim::evaluateLegendrePolynomial(const int nbasis, const int nt, const double* coeff, Vector& params)
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

double uniformRand(double min, double max)
{
//  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(2.0); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> uniform_rand(0, 1);
  return min + (max - min)*uniform_rand(gen);
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

Matrix getDCTMatrix(int N)
{
  const double scale_row0 = 1.0 / std::sqrt(N);
  const double scale_rowk = std::sqrt(2.0) / std::sqrt(N);

  Matrix A(N, N, 0.0);

  //A(0,:)
  for (int n = 0; n < N; ++n)
    A(0,n) = scale_row0; //row k = 0

  //A(1:N,:)
  for (int k = 1; k < N; ++k)
    for (int n = 0; n < N; ++n)
      A(k,n) = std::cos(M_PI*(n + 0.5)*k / N) * scale_rowk;

  return A;
}

std::string
getNLOPTResultDescription(nlopt::result resultcode)
{
  switch (resultcode)
  {
    case nlopt::result::FAILURE:
      return "Failed - Generic result code";
      break;

    case nlopt::result::INVALID_ARGS:
      return "Failed - Invalid arguments";
      break;

    case nlopt::result::OUT_OF_MEMORY:
      return "Failed - Out of memory";
      break;

    case nlopt::result::ROUNDOFF_LIMITED:
      return "Failed - Round-off limited";
      break;

    case nlopt::result::FORCED_STOP:
      return "Failed - Forcefully stopped";
      break;

    case nlopt::result::SUCCESS:
      return "Success - Generic result code";
      break;

    case nlopt::result::STOPVAL_REACHED:
      return "Success - Stop value reached";
      break;

    case nlopt::result::FTOL_REACHED:
      return "Success - Relative f-tolerance reached";
      break;

    case nlopt::result::XTOL_REACHED:
      return "Success - Relative x-tolerance reached";
      break;

    case nlopt::result::MAXEVAL_REACHED:
      return "Success - Maximum evaluation count reached";
      break;

    case nlopt::result::MAXTIME_REACHED:
      return "Success - Maximum time reached";
      break;

    default:
      break;
  }

  return "Unknown result code";
}
