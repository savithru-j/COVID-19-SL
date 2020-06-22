#include <iostream>
#include <iomanip>
#include <cassert>

#include "OptimizerGlobalBasis.h"
#include "Simulator.h"

OptimizerGlobalBasis::OptimizerGlobalBasis(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                                 const Vector& quarantine_input, int num_basis_,
                                 double wconf, double wrecov, double wfatal,
                                 int max_iter_per_pass_, int max_passes_, int seed) :
    pop_observed(pop_observed_), pop_init(pop_init_),
    nt_opt(pop_observed.getNumDays()), num_basis(num_basis_),
    weight_conf(wconf), weight_recov(wrecov), weight_fatal(wfatal),
    max_iter_per_pass(max_iter_per_pass_), max_passes(max_passes_),
    params(pop_observed.getNumDays(), 0, quarantine_input),
    rand_engine(seed), uniform_rand(0,1)
{
  if (pop_observed.N != pop_init.N)
    throwError("Initial population mismatch!");

  transform_matrix = getDCTMatrix(params.nt_hist, num_basis);
  inv_transform_matrix = getInverseDCTMatrix(params.nt_hist, num_basis);

  param_bounds = getParameterBounds(nt_opt, num_basis);
  param_vec.resize(param_bounds.size());
  copyParam2Vector(params, param_vec);
}


void
OptimizerGlobalBasis::optimizeParametersNLOPT()
{
  std::cout << "Optimizing parameters for " << params.nt_hist << " days..." << std::endl;
  std::cout << "No. of parameters: " << nDim() << std::endl;

  f_eval_count = 0;

  //nlopt optimizer object
  nlopt::opt opt(nlopt::LN_COBYLA, nDim());

  opt.set_min_objective( getCostNLOPT, reinterpret_cast<void*>(this) );

  const int num_constraints = 2*4; //2*4*(num_basis+1); //(min and max)*(4 params)*(at (num_basis+1) locations)
  std::cout << "No. of constraints: " << num_constraints << std::endl;

#if 0 //Numerically check constraint gradient
  randomizeParameters();

  std::vector<double> con0(num_constraints, 0.0);
  std::vector<double> con1(num_constraints, 0.0);
  std::vector<double> grad0(num_constraints*nDim(), 0.0);
  std::vector<double> gradFD(num_constraints*nDim(), 0.0);
  std::vector<double> grad_dummy(num_constraints*nDim(), 0.0);
  getConstraintsDCT(con0.data(), grad0.data());

  for (int j = 0; j < nDim(); ++j)
  {
    param_vec[j] += 1.0;
    getConstraintsDCT(con1.data(), grad_dummy.data());

    for (int i = 0; i < num_constraints; ++i)
      gradFD[i*nDim() + j] = con1[i] - con0[i];

    param_vec[j] -= 1.0;
  }

  for (int i = 0; i < con0.size(); ++i)
    std::cout << i << ": " << con0[i] << std::endl;

  for (int i = 0; i < num_constraints; ++i)
  {
    for (int j = 0; j < nDim(); ++j)
    {
      int ind = i*nDim() + j;
//      if (std::abs(grad0[ind] - gradFD[ind]) > 1e-13)
        std::cout << i << ", " << j << ": " << grad0[ind] << ", " << gradFD[ind] << ", " << grad0[ind] - gradFD[ind] << std::endl;
    }
  }
  throw;
#endif

  std::vector<double> constraint_tol(num_constraints, 0.0);
  opt.add_inequality_mconstraint(getConstraintsNLOPT, reinterpret_cast<void*>(this), constraint_tol);

  //stop when every parameter changes by less than the tolerance multiplied by the absolute value of the parameter.
  opt.set_xtol_rel(1e-6);

  //stop when the objective function value changes by less than the tolerance multiplied by the absolute value of the function value
  opt.set_ftol_rel(1e-8);

  //stop when the maximum number of function evaluations is reached
  opt.set_maxeval(max_iter_per_pass);

  int num_stored = 0;

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

    double penalty = getPenalty();
    std::cout << "Penalty: " << penalty << std::endl;

    if (cost_opt < cost_min.back())
    {
      updateOptimalSolution(cost_opt, subcosts_dummy, x);
      num_stored++;
    }

//    param_vec = x;
//    Vector con(num_constraints);
//    getConstraintsDCT(con.getDataVector().data(), nullptr);
//    std::cout << "Con: " << con << std::endl;
//    std::cout << "Params: " << param_vec << std::endl;

    randomizeParameters();
    x = param_vec.getDataVector();
  }


  for (int i = num_stored; i < NUM_RESULTS; ++i)
  {
    optimal_param_vec[i] = optimal_param_vec[num_stored-1];
    cost_min[i] = cost_min[num_stored-1];
  }

  std::cout << std::endl;
  std::cout << "Minimum costs (relative): " << cost_min << std::endl;
  std::cout << "Cost function evaluation count: " << f_eval_count << std::endl;
}

double
OptimizerGlobalBasis::getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
  OptimizerGlobalBasis* opt = reinterpret_cast<OptimizerGlobalBasis*>(data);

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

void
OptimizerGlobalBasis::getConstraintsNLOPT(
    unsigned m, double *constraints, unsigned n, const double* x, double* grad, void* data)
{
  OptimizerGlobalBasis* opt = reinterpret_cast<OptimizerGlobalBasis*>(data);

  assert((int) m == 2*4); //*(opt->num_basis+1));
  assert((int) n == opt->nDim());

  for (unsigned int i = 0; i < n; ++i)
    opt->param_vec[i] = x[i];
  opt->copyVector2Param(opt->param_vec, opt->params);

  std::fill(constraints, constraints + m, 0.0);

  if (grad)
    std::fill(grad, grad + m*n, 0.0);

//  opt->getConstraints(constraints, grad);
  opt->getConstraintsDCT(constraints, grad);
}

std::pair<double, std::array<double, 3>>
OptimizerGlobalBasis::getCost()
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

  const double cost = sub_costs[0] + sub_costs[1] + sub_costs[2]; // + getPenalty();

  f_eval_count++;
  return {cost, sub_costs};
}

double
OptimizerGlobalBasis::getCostGradient(std::vector<double>& grad)
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
OptimizerGlobalBasis::getConstraints(double* constraints, double* grad)
{
  const int num_cpts = num_basis+1;  //no. of points where constraints are enforced
  const int off = 4*num_cpts; //index offset for maximum constraints
  const int ndim = nDim();

  Vector phi(num_basis, 0.0);

  /* Constraint equation form:
   *    sum( u[i]*phi[i] ) - u_max <= 0.0
   *   -sum( u[i]*phi[i] ) + u_min <= 0.0
   */

  for (int i = 0; i < num_cpts; ++i)
  {
    const double x = 2.0*(i / (double)(num_cpts - 1)) - 1.0; //in range [-1, 1]
    getBasisLegendre(x, phi);

    constraints[             i] = param_bounds[          0].min;
    constraints[  num_cpts + i] = param_bounds[  num_basis].min;
    constraints[2*num_cpts + i] = param_bounds[2*num_basis].min;
    constraints[3*num_cpts + i] = param_bounds[3*num_basis].min;

    constraints[off +              i] = -param_bounds[          0].max;
    constraints[off +   num_cpts + i] = -param_bounds[  num_basis].max;
    constraints[off + 2*num_cpts + i] = -param_bounds[2*num_basis].max;
    constraints[off + 3*num_cpts + i] = -param_bounds[3*num_basis].max;


    for (int j = 0; j < num_basis; ++j)
    {
      constraints[             i] -= param_vec[              j] * phi[j];
      constraints[  num_cpts + i] -= param_vec[  num_basis + j] * phi[j];
      constraints[2*num_cpts + i] -= param_vec[2*num_basis + j] * phi[j];
      constraints[3*num_cpts + i] -= param_vec[3*num_basis + j] * phi[j];

      constraints[off +              i] += param_vec[              j] * phi[j];
      constraints[off +   num_cpts + i] += param_vec[  num_basis + j] * phi[j];
      constraints[off + 2*num_cpts + i] += param_vec[2*num_basis + j] * phi[j];
      constraints[off + 3*num_cpts + i] += param_vec[3*num_basis + j] * phi[j];

      grad[(             i)*ndim +               j] = -phi[j];
      grad[(  num_cpts + i)*ndim +   num_basis + j] = -phi[j];
      grad[(2*num_cpts + i)*ndim + 2*num_basis + j] = -phi[j];
      grad[(3*num_cpts + i)*ndim + 3*num_basis + j] = -phi[j];

      grad[(off +              i)*ndim +               j] = phi[j];
      grad[(off +   num_cpts + i)*ndim +   num_basis + j] = phi[j];
      grad[(off + 2*num_cpts + i)*ndim + 2*num_basis + j] = phi[j];
      grad[(off + 3*num_cpts + i)*ndim + 3*num_basis + j] = phi[j];
    }
  }
}

void
OptimizerGlobalBasis::getConstraintsDCT(double* constraints, double* grad)
{
  const int ndim = nDim();
  const int nt = params.nt_hist;

  const double f0 = std::sqrt(1.0/nt);
  const double f1 = std::sqrt(2.0/nt);

  std::array<double,4> xmin = {0, 0, 0, 1};
  std::array<double,4> xmax = {2, 1, 1, 1};

  for (int p = 0; p < 4; ++p)
  {
    constraints[2*p    ] =  xmin[p] - f0 * param_vec[p*num_basis];
    constraints[2*p + 1] = -xmax[p] + f0 * param_vec[p*num_basis];

    for (int j = 1; j < num_basis; ++j)
    {
      constraints[2*p    ] += f1 * std::abs(param_vec[p*num_basis + j]);
      constraints[2*p + 1] += f1 * std::abs(param_vec[p*num_basis + j]);
    }
  }

  if (grad)
  {
    for (int p = 0; p < 4; ++p)
    {
      grad[(2*p    )*ndim + p*num_basis] = -f0;
      grad[(2*p + 1)*ndim + p*num_basis] =  f0;

      for (int j = 1; j < num_basis; ++j)
      {
        grad[(2*p    )*ndim + p*num_basis + j] = (param_vec[p*num_basis + j] >= 0) ? f1 : -f1;
        grad[(2*p + 1)*ndim + p*num_basis + j] = (param_vec[p*num_basis + j] >= 0) ? f1 : -f1;
      }
    }
  }

}

double
OptimizerGlobalBasis::getPenalty()
{
  std::array<double,4> xmin = {0, 0, 0, 1};
  std::array<double,4> xmax = {2, 0, 0, 1};

  double cost = 0.0;
  for (int i = 0; i < params.nt_hist; ++i)
  {
    cost += std::max(0.0, xmin[0] - params.betaN[i]) + std::max(0.0, params.betaN[i] - xmax[0]);
    cost += std::max(0.0, xmin[1] - params.c0[i]) + std::max(0.0, params.c0[i] - xmax[1]);
    cost += std::max(0.0, xmin[2] - params.c1[i]) + std::max(0.0, params.c1[i] - xmax[2]);
    cost += std::max(0.0, xmin[3] - params.c2[i]) + std::max(0.0, params.c2[i] - xmax[3]);
  }
  return cost;
}

void
OptimizerGlobalBasis::updateOptimalSolution(
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

std::vector<ParamBound>
OptimizerGlobalBasis::getParameterBounds(int nt, int nbasis)
{
  const int m = 4*nbasis + 11;
  std::vector<ParamBound> bounds(m);

  const double delta = 1e-4;
//  const double a = 1e2;
  const double f0 = std::sqrt(nt);
  const double f1 = std::sqrt(2.0*nt);

  const double beta_max = 2.0;
  //Bounds for "average" basis coefficients
  bounds[       0] = ParamBound(0*f0, beta_max*f0, delta); //betaN
  bounds[  nbasis] = ParamBound(0*f0, 1*f0, delta); //c0 = ce
  bounds[2*nbasis] = ParamBound(0*f0, 1*f0, delta); //c1
  bounds[3*nbasis] = ParamBound(1*f0, 1*f0, delta); //c2 = c3

  for (int i = 1; i < nbasis; ++i)
  {
    bounds[           i] = ParamBound(-beta_max*f1, beta_max*f1, delta); //betaN
    bounds[  nbasis + i] = ParamBound(-f1, f1, delta); //c0 = ce
    bounds[2*nbasis + i] = ParamBound(-f1, f1, delta); //c1
    bounds[3*nbasis + i] = ParamBound(-f1*0, f1*0, delta); //c2 = c3
  }


  const int off = 4*nbasis;
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
OptimizerGlobalBasis::copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int m = 4*num_basis + 11;
  if (v.m() != m)
    throwError("copyParam2Vector - inconsistent dimensions!");

  std::fill(v.begin(), v.end(), 0.0);

  for (int i = 0; i < num_basis; ++i)
  {
    for (int j = 0; j < params.nt_hist; ++j)
    {
      v[              i] += transform_matrix(i,j) * params.betaN[j];
      v[  num_basis + i] += transform_matrix(i,j) * params.c0[j];
      v[2*num_basis + i] += transform_matrix(i,j) * params.c1[j];
      v[3*num_basis + i] += transform_matrix(i,j) * params.c2[j];
    }
  }

  const int off = 4*num_basis;
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
OptimizerGlobalBasis::copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 4*num_basis + 11;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < params.nt_hist; ++i)
  {
    params.betaN[i] = 0.0;
    params.c0[i] = 0.0;
    params.c1[i] = 0.0;
    params.c2[i] = 0.0;

    for (int j = 0; j < num_basis; ++j)
    {
      params.betaN[i] += inv_transform_matrix(i,j) * v[              j];
      params.c0[i]    += inv_transform_matrix(i,j) * v[  num_basis + j];
      params.c1[i]    += inv_transform_matrix(i,j) * v[2*num_basis + j];
      params.c2[i]    += inv_transform_matrix(i,j) * v[3*num_basis + j];
    }
  }

#if 0
  Vector phi(num_basis, 0.0);

  for (int i = 0; i < nt; ++i)
  {
    const double x = 2.0*(i / (double)(nt - 1)) - 1.0; //in range [-1, 1]
    getBasisLegendre(x, phi);

    params.betaN[i] = 0.0;
    params.c0[i] = 0.0;
    params.c1[i] = 0.0;
    params.c2[i] = 0.0;

    for (int j = 0; j < num_basis; ++j)
    {
      params.betaN[i] += v[              j] * phi[j];
      params.c0[i]    += v[  num_basis + j] * phi[j];
      params.c1[i]    += v[2*num_basis + j] * phi[j];
      params.c2[i]    += v[3*num_basis + j] * phi[j];
    }
  }
#endif

  params.ce = params.c0;
  params.c3 = params.c2;

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

  const int off = 4*num_basis;
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
OptimizerGlobalBasis::getBasisLegendre(const double& x, Vector& phi)
{
  switch (phi.size())
  {
  case 1:
    phi[0] = 1.0;
    break;

  case 2:
    phi[0] = 1.0;
    phi[1] = x;
    break;

  case 3:
    phi[0] = 1.0;
    phi[1] = x;
    phi[2] = 1.5*x*x - 0.5;
    break;

  case 4:
    phi[0] = 1.0;
    phi[1] = x;
    phi[2] = 1.5*x*x - 0.5;
    phi[3] = 0.5*x*(5.0*x*x - 3.0);
    break;

  default:
    throwError("getBasisLegendre - Unsupported basis number!");
  }
}
