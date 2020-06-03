#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <nlopt.hpp>

#include "Population.h"
#include "LinearAlgebra.h"

struct ParamBound {

  ParamBound() = default;
  ParamBound(double min_, double max_, double step_) : min(min_), max(max_), step(step_) {}

  double min = 0.0, max = 1.0, step = 1e-4;
};

double uniformRand(double min = 0.0, double max = 1.0);

Matrix getHaarMatrix(int m);
Matrix getDCTMatrix(int N);

std::string getNLOPTResultDescription(nlopt::result resultcode);

struct Optimizer
{
  static constexpr int NUM_RESULTS = 5;  //no. of optimal results to store (best to worst)

  Optimizer(const ObservedPopulation& pop_observed_, const Population& pop_init_,
            double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1, double wreg_ = 0.01,
            int max_iter_per_pass_ = 1000, int max_passes_ = 1);

  inline int nDim() const { return param_vec.size(); };

  inline void randomizeParameters(const double energy = 1.0)
  {
    for (int i = 0; i < param_vec.m(); ++i)
    {
      const double delta = energy * (param_bounds[i].max - param_bounds[i].min) * uniformRand(-1, 1);
      param_vec[i] = std::min(std::max(param_vec[i] + delta, param_bounds[i].min), param_bounds[i].max);
    }
    copyVector2Param(param_vec, params);
  }

  void optimizeParameters();
  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population& pop_init;
  double weight_conf, weight_recov, weight_fatal, weight_reg;
  int max_iter_per_pass = 1000;
  int max_passes = 1;

  const int t_buffer = 0; //no. of non-optimized days at end
  const int nt_opt; //number of days to optimize parameters for

  double cost_reduction_tol = 1e-4;
  double min_eta = 1e-5;

  std::vector<ParamBound> param_bounds;
  Vector param_vec; //current solution vector
  ModelParams params;

  Matrix reg_matrix = Matrix(0,0);

  int f_eval_count = 0;
  int nlopt_iter = 0;

  std::array<Vector,NUM_RESULTS> optimal_param_vec;
  std::array<double,NUM_RESULTS> cost_rel_min;
  std::array<std::array<double,6>,NUM_RESULTS> sub_costs_min;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  static void copyParam2Vector(const ModelParams& params, Vector& v);
  static void copyVector2Param(const Vector& v, ModelParams& params);

protected:

  std::pair<double, std::array<double, 6>> getCost();

  double getCostGradient(std::vector<double>& grad);
  double getCostGradient(Vector& grad) { return getCostGradient(grad.getDataVector()); }

  static std::vector<ParamBound> getParameterBounds(int nt);

  void updateOptimalSolution(const double& cost_rel, const std::array<double,6>& sub_costs,
                             const Vector& param_vec);

  double limitUpdate(Vector& dparam_vec);
};


struct OptimizerLowDim
{
  static constexpr int NUM_RESULTS = 5;  //no. of optimal results to store (best to worst)

  OptimizerLowDim(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                  int num_basis_, double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1,
                  int max_iter_per_pass_ = 1000, int max_passes_ = 1);

  inline int nDim() const { return param_vec.size(); };

#if 0
  inline void randomizeParameters()
  {
    for (int i = 0; i < 5; ++i)
    {
      param_vec[i*num_basis] = uniformRand(param_bounds[i*num_basis].min, param_bounds[i*num_basis].max);
      for (int j = 1; j < num_basis; ++j)
        param_vec[i*num_basis + j] = 0.0;
    }

    const int off = 5*num_basis;
    for (int i = off; i < param_vec.m(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);

    copyVector2Param(param_vec, params);
  }
#else
  inline void randomizeParameters()
  {
    const int num_nodes = (int)(nt_opt/num_basis) + 1;
    for (int i = 0; i < 5; ++i)
    {
      param_vec[i*num_nodes] = uniformRand(param_bounds[i*num_nodes].min, param_bounds[i*num_nodes].max);
      for (int j = 1; j < num_nodes; ++j)
        param_vec[i*num_nodes + j] = param_vec[i*num_nodes];
    }

    const int off = 5*num_nodes;
    for (int i = off; i < param_vec.m(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);

    copyVector2Param(param_vec, params);
  }
#endif

  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population& pop_init;
  const int nt_opt, num_basis;
  double weight_conf, weight_recov, weight_fatal;
  int max_iter_per_pass = 1000;
  int max_passes = 1;

  std::vector<ParamBound> param_bounds;
  Vector param_vec; //current solution vector
  ModelParams params;

  int f_eval_count = 0;
  int nlopt_iter = 0;

  std::array<Vector,NUM_RESULTS> optimal_param_vec;
  std::array<double,NUM_RESULTS> cost_min;
  std::array<std::array<double,3>,NUM_RESULTS> sub_costs_min;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  void copyParam2Vector(const ModelParams& params, Vector& v);
  void copyVector2Param(const Vector& v, ModelParams& params);

protected:

  std::pair<double, std::array<double, 3>> getCost();

  double getCostGradient(std::vector<double>& grad);
  double getCostGradient(Vector& grad) { return getCostGradient(grad.getDataVector()); }

  static std::vector<ParamBound> getParameterBounds(int nt, int num_basis);

  static void evaluateLegendrePolynomial(const int nbasis, const int nt, const double* coeff, Vector& params);

  void updateOptimalSolution(const double& cost_rel, const std::array<double,3>& sub_costs,
                             const Vector& param_vec);
};

#endif
