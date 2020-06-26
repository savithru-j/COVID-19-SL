#include <iostream>
#include <iomanip>
#include <limits>

#include "OptimizerFull.h"
#include "Simulator.h"

OptimizerFull::OptimizerFull(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                             const Vector& quarantine_input,
                             double wconf, double wrecov, double wfatal, double wreg,
                             int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal), weight_reg(wreg),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    nt_opt(pop_observed.getNumDays() - t_buffer),
    params(nt_opt, t_buffer, quarantine_input),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  if (weight_reg != 0.0)
  {
    reg_matrix = getHaarMatrix(nt_opt);
//    reg_matrix = getDCTMatrix(nt_opt, nt_opt);
  }

  param_bounds = getParameterBounds(nt_opt);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}

void
OptimizerFull::optimizeParameters()
{
  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;

  f_eval_count = 0;

  double jump_energy = 1.0;

  Vector param_vec0 = param_vec;
  optimal_param_vec.fill(param_vec);

  double cost0 = -1;
  cost_min.fill(1.0); //relative cost

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

    if (cost_rel < cost_min.back())
    {
      updateOptimalSolution(cost_rel, cost.second, param_vec);
    }

    jump_energy = 1.0 - pass / (double) (max_passes - 1.0);
    std::cout << "Jump energy: " << jump_energy << std::endl;

    param_vec = optimal_param_vec[0]; //goto current optimal solution
    randomizeParameters(jump_energy);

  } //passes

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Minimum sub-costs [sol-error, beta-reg, c0-reg, c1-reg]:" << std::endl;
  for (std::size_t i = 0; i < NUM_RESULTS; i++)
      std::cout << "  " << sub_costs_min[i] << std::endl;

  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;

//  for (std::size_t i = 0; i < NUM_RESULTS; i++)
//    std::cout << "Optimal params[" << i << "]: " << param_vec_opt[i] << std::endl << std::endl;
}

void
OptimizerFull::optimizeParametersNLOPT()
{
  std::cout << "Optimizing parameters for " << nt_opt << " days..." << std::endl;
  std::cout << "No. of parameters: " << nDim() << std::endl;

  f_eval_count = 0;

  //nlopt optimizer object
  nlopt::opt opt(nlopt::LD_LBFGS, nDim());
//  nlopt::opt opt(nlopt::LN_COBYLA, nDim());

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
  std::array<double,6> subcosts_dummy;

  std::cout << std::scientific << std::setprecision(4);

  for (int pass = 0; pass < max_passes; ++pass)
  {
    std::cout << std::endl << "Pass " << (pass+1) << "/" << max_passes << " ----------------" << std::endl;

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

    randomizeParameters(1.0);
    x = param_vec.getDataVector();
  }

  for (int i = max_passes; i < NUM_RESULTS; ++i)
  {
    optimal_param_vec[i] = optimal_param_vec[max_passes-1];
    cost_min[i] = cost_min[max_passes-1];
  }

  std::cout << std::endl;
  std::cout << "Minimum costs: " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerFull::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerFull* opt = reinterpret_cast<OptimizerFull*>(data);

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
OptimizerFull::getCost()
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
OptimizerFull::getCostGradient(std::vector<double>& grad)
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
OptimizerFull::updateOptimalSolution(
    const double& cost, const std::array<double,6>& sub_costs, const Vector& param_vec_cur)
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

double
OptimizerFull::limitUpdate(Vector& dparam_vec)
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
OptimizerFull::copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int nt = params.nt_hist;
  const int m = 2*nt;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    v[     i] = params.betaN[i];
    v[nt + i] = params.c1[i];
  }
//  const int off = 5*nt;
//  v[off  ] = params.T_incub0;
//  v[off+1] = params.T_incub1;
//  v[off+2] = params.T_asympt;
//  v[off+3] = params.T_mild;
//  v[off+4] = params.T_severe;
//  v[off+5] = params.T_icu;
//  v[off+6] = params.f;
//  v[off+7] = params.frac_recover_I1;
//  v[off+8] = params.frac_recover_I2;
//  v[off+9] = params.CFR;
//  v[off+10] = params.T_discharge;
}

void
OptimizerFull::copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 2*nt;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[     i];
    params.c1[i]    = v[nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = params.betaN[nt-1];
    params.c1[i]    = params.c1[nt-1];
  }

//  const int off = 5*nt;
//  params.T_incub0        = v[off  ];
//  params.T_incub1        = v[off+1];
//  params.T_asympt        = v[off+2];
//  params.T_mild          = v[off+3];
//  params.T_severe        = v[off+4];
//  params.T_icu           = v[off+5];
//  params.f               = v[off+6];
//  params.frac_recover_I1 = v[off+7];
//  params.frac_recover_I2 = v[off+8];
//  params.CFR             = v[off+9];
//  params.T_discharge     = v[off+10];
}

std::vector<ParamBound>
OptimizerFull::getParameterBounds(int nt)
{
  const int m = 2*nt;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;

  for (int i = 0; i < nt; ++i)
  {
    bounds[     i] = ParamBound(0, 2, delta); //betaN
    bounds[nt + i] = ParamBound(0, 1, delta); //c1
  }
//  const int off = 5*nt;
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
  return bounds;
}
