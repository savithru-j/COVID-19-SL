#ifndef OPTIMIZERPIECEWISELINEAR_H
#define OPTIMIZERPIECEWISELINEAR_H

#include <random>

#include "Optimization.h"
#include "Population.h"
#include "LinearAlgebra.h"

struct OptimizerPiecewiseLinear
{
  static constexpr int NUM_RESULTS = 40;  //no. of optimal results to store (best to worst)

  OptimizerPiecewiseLinear(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                           const Vector& quarantine_input, int num_basis_,
                           double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1,
                           int max_iter_per_pass_ = 1000, int max_passes_ = 1, int seed = 1);

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
#elif 1
  inline void randomizeParameters()
  {
    const int num_nodes = (int)(nt_opt/num_basis) + 1;
    for (int i = 0; i < 4; ++i)
    {
//       param_vec[i*num_nodes] = uniformRand(param_bounds[i*num_nodes].min, param_bounds[i*num_nodes].max);
//       for (int j = 1; j < num_nodes; ++j)
//         param_vec[i*num_nodes + j] = param_vec[i*num_nodes];

      for (int j = 0; j < num_nodes; ++j)
        param_vec[i*num_nodes + j] = uniformRand(param_bounds[i*num_nodes+j].min, param_bounds[i*num_nodes+j].max);

    }

    const int off = 4*num_nodes;
    for (int i = off; i < param_vec.m(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);

    copyVector2Param(param_vec, params);
  }
#else
  inline void randomizeParameters()
  {
    const int num_nodes = (int)(nt_opt/num_basis) + 1;
    param_vec[0] = uniformRand(param_bounds[0].min, param_bounds[0].max);
    for (int j = 1; j < num_nodes; ++j)
      param_vec[j] = param_vec[0];

    const int off = num_nodes;
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

  std::mt19937 rand_engine; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> uniform_rand;

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

  inline double uniformRand(double min = 0.0, double max = 1.0)
  {
    return min + (max - min)*uniform_rand(rand_engine);
  }
};

#endif
