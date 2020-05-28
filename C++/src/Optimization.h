#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "Population.h"
#include "LinearAlgebra.h"

struct ParamBound {

  ParamBound() = default;
  ParamBound(double min_, double max_, double step_) : min(min_), max(max_), step(step_) {}

  double min = 0.0, max = 1.0, step = 1e-4;
};

void copyParam2Vector(const ModelParams& params, Vector& v);
void copyVector2Param(const Vector& v, ModelParams& params);

std::vector<ParamBound> getParameterBounds(int nt);

double uniformRand(double min = 0.0, double max = 1.0);

Matrix getHaarMatrix(int m);



struct Optimizer
{
  static constexpr int NUM_RESULTS = 5;  //no. of optimal results to store (best to worst)

  Optimizer(const ObservedPopulation& pop_observed_, const Population& pop_init_, double reg_weight_ = 0.1);

//  inline int nOptDays() const { return nt_opt; }
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

  int t_buffer = 0*5; //no. of non-optimized days at end

  double cost_reduction_tol = 1e-4;
  double min_eta = 1e-5;
  int max_iter_per_pass = 200;
  int max_passes = 10;

  const ObservedPopulation& pop_observed;
  const Population& pop_init;

  double reg_weight;
  int nt_opt; //number of days to optimize parameters for

  std::vector<ParamBound> param_bounds;
  Vector param_vec; //current solution vector
  ModelParams params;

  Matrix reg_matrix = Matrix(0,0);

  int f_eval_count = 0;

  std::array<Vector,NUM_RESULTS> optimal_param_vec;
  std::array<double,NUM_RESULTS> cost_rel_min;
  std::array<std::array<double,4>,NUM_RESULTS> sub_costs_min;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

protected:

  std::pair<double, std::array<double, 4>> getCost();

  void getCostGradient(std::vector<double>& grad);
  void getCostGradient(Vector& grad) { getCostGradient(grad.getDataVector()); }

  void updateOptimalSolution(const double& cost_rel, const std::array<double,4>& sub_costs,
                             const Vector& param_vec);

  double limitUpdate(Vector& dparam_vec);
};

#endif
