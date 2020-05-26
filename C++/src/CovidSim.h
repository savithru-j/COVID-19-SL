#ifndef COVIDSIM_H
#define COVIDSIM_H

#include <random>

#include "Population.h"
#include "LinearAlgebra.h"

struct ParamBound {

  ParamBound() = default;
  ParamBound(double min_, double max_, double step_) : min(min_), max(max_), step(step_) {}

  double min = 0.0, max = 1.0, step = 1e-4;
};

struct OptimizationInfo
{
  static constexpr int NUM_RESULTS = 5;  //no. of optimal results to store (best to worst)
  static constexpr int T_END_BUFFER = 5; //no. of non-optimized days at end

  static constexpr double COST_REDUCTION_TOL = 1e-4;
  static constexpr double MIN_ETA = 1e-5;
  static constexpr int MAX_ITER_PER_PASS = 50;
  static constexpr int MAX_PASSES = 100;

  Matrix reg_matrix = Matrix(0,0);
  bool regularize = true;
  double cost_error_scaling = 1.0;
  int f_eval_count = 0;

  std::array<Vector,NUM_RESULTS> param_vec_opt;
  std::array<double,NUM_RESULTS> cost_rel_min;
  std::array<double,NUM_RESULTS> cost_raw_min;
  std::array<std::array<double,4>,NUM_RESULTS> sub_costs_min;

  void updateOptimalSolution(const double& cost_rel, const double& cost_raw,
                             const std::array<double,4>& sub_costs, const Vector& param_vec)
  {
    int min_ind = 0;
    for (min_ind = 0; min_ind < NUM_RESULTS; ++min_ind)
      if (cost_rel < cost_rel_min[min_ind])
        break;

    for (int i = NUM_RESULTS-1; i > min_ind; i--)
    {
      cost_rel_min[i] = cost_rel_min[i-1];
      cost_raw_min[i] = cost_raw_min[i-1];
      sub_costs_min[i] = sub_costs_min[i-1];
      param_vec_opt[i] = param_vec_opt[i-1];
    }
    cost_rel_min[min_ind] = cost_rel;
    cost_raw_min[min_ind] = cost_raw;
    sub_costs_min[min_ind] = sub_costs;
    param_vec_opt[min_ind] = param_vec;
  }
};


int main();

std::vector<Population> predictModel(const ModelParams& params, const Population& pop_init);

OptimizationInfo optimizeParameters(const ObservedPopulation& pop_observed, const Population& pop_init);

std::pair<double, std::array<double, 4>>
getCost(const ModelParams& params, const Population& pop_init,
        const ObservedPopulation& pop_observed, OptimizationInfo& optinfo);

Vector getCostGradient(const ModelParams& params, const Population& pop_init,
                       const ObservedPopulation& pop_observed, const std::vector<ParamBound>& bounds,
                       OptimizationInfo& optinfo);

void copyParam2Vector(const ModelParams& params, Vector& v);
void copyVector2Param(const Vector& v, ModelParams& params);

std::vector<ParamBound> getParameterBounds(int nt);

double uniformRand(double min = 0.0, double max = 1.0);

void randomizeParameterVector(const std::vector<ParamBound>& bounds, const double energy, Vector& param_vec);

double limitUpdate(const std::vector<ParamBound>& bounds, const Vector& param_vec, Vector& dparam_vec);

Matrix getHaarMatrix(int m);

#endif
