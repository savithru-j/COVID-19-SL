#ifndef OPTIMIZERFUNCTIONAL_H
#define OPTIMIZERFUNCTIONAL_H

#include <random>

#include "Optimization.h"
#include "Population.h"
#include "LinearAlgebra.h"

struct OptimizerFunctional
{
  OptimizerFunctional(const ObservedPopulation& pop_observed_, const Population<double>& pop_init_,
                      const Vector<double>& quarantine_input, const int T_avg_, const int T_unopt_,
                      int max_iter_per_pass_ = 1000, int max_passes_ = 1, int seed = 1);

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

  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population<double>& pop_init;
  int T_avg, T_unopt;
  int max_iter_per_pass;
  int max_passes;

  double cost_reduction_tol = 1e-4;
//  double min_eta = 1e-5;

  std::vector<ParamBound> param_bounds;
  Vector<double> param_vec; //current solution vector
  ModelParams<double> params;

  int f_eval_count = 0;
  int nlopt_iter = 0;

  int num_results; //no. of optimal results to store (best to worst)
  std::vector<Vector<double>> optimal_param_vec;
  std::vector<double> cost_min;
  std::vector<std::array<double,6>> sub_costs_min;

  std::mt19937 rand_engine; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> uniform_rand;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  static void copyParam2Vector(const ModelParams<double>& params, Vector<double>& v);
  static void copyVector2Param(const Vector<double>& v, ModelParams<double>& params);

protected:

  std::pair<double, std::array<double, 6>> getCost();

  double getCostGradient(std::vector<double>& grad);
  double getCostGradient(Vector<double>& grad) { return getCostGradient(grad.getDataVector()); }

  static std::vector<ParamBound> getParameterBounds(int nt);

  void updateOptimalSolution(const double& cost, const std::array<double,6>& sub_costs,
                             const Vector<double>& param_vec);

  double limitUpdate(Vector<double>& dparam_vec);

  inline double uniformRand(double min = 0.0, double max = 1.0)
  {
    return min + (max - min)*uniform_rand(rand_engine);
  }
};

#endif
