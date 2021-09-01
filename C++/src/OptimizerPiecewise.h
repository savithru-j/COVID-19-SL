#ifndef OPTIMIZERPIECEWISE_H
#define OPTIMIZERPIECEWISE_H

#include <random>

#include "Optimization.h"
#include "Population.h"
#include "LinearAlgebra.h"

struct OptimizerPiecewise
{
//  static constexpr int NUM_RESULTS = 40;  //no. of optimal results to store (best to worst)
  static constexpr bool OPTIMIZE_C0 = true;
  static constexpr bool OPTIMIZE_C1 = true;
  static constexpr bool OPTIMIZE_C2 = false;
  static constexpr std::array<int,3> IFR_SEG_STARTS = {0, 182};//, 273};

  OptimizerPiecewise(const ObservedPopulation& pop_observed_, const Population& pop_init_,
                           const Vector& quarantine_input, int interval_size_, bool linear_basis_ = false,
                           double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1,
                           int max_iter_per_pass_ = 1000, int max_passes_ = 1, int seed = 1);

  inline int nDim() const { return param_vec.size(); };

  inline void randomizeParameters()
  {
    const int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
    constexpr int num_c_params = OPTIMIZE_C0 + OPTIMIZE_C1 + OPTIMIZE_C2;
    for (int i = 0; i < (1+num_c_params); ++i)
    {
#if 1 //Constant random solutions
      param_vec[i*num_nodes] = uniformRand(param_bounds[i*num_nodes].min, param_bounds[i*num_nodes].max);
      for (int j = 1; j < num_nodes; ++j)
        param_vec[i*num_nodes + j] = param_vec[i*num_nodes];
#else
      //Piecewise constant/linear random solutions
      for (int j = 0; j < num_nodes; ++j)
        param_vec[i*num_nodes + j] = uniformRand(param_bounds[i*num_nodes+j].min, param_bounds[i*num_nodes+j].max);
#endif
    }

    //Random constant solution for IFR (same value for all segments)
    int off = (1+num_c_params)*num_nodes;
    param_vec[off] = uniformRand(param_bounds[off].min, param_bounds[off].max);
    for (std::size_t i = 1; i < IFR_SEG_STARTS.size(); ++i)
      param_vec[off + i] = param_vec[off];

    off += IFR_SEG_STARTS.size();
    for (int i = off; i < param_vec.m(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);

    copyVector2Param(param_vec, params);
  }

  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population& pop_init;
  const int nt_opt, interval_size;
  const bool linear_basis = false;
  double weight_conf, weight_recov, weight_fatal;
  int max_iter_per_pass = 1000;
  int max_passes = 1;

  std::vector<ParamBound> param_bounds;
  Vector param_vec; //current solution vector
  ModelParams params;

  int f_eval_count = 0;
  int nlopt_iter = 0;

  std::vector<Vector> optimal_param_vec;
  std::vector<double> cost_min;
  std::vector<std::array<double,3>> sub_costs_min;

  std::mt19937 rand_engine; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> uniform_rand;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  void copyParam2Vector(const ModelParams& params, Vector& v);
  void copyVector2Param(const Vector& v, ModelParams& params);

protected:

  std::pair<double, std::array<double, 3>> getCost();

  double getCostGradient(std::vector<double>& grad);
  double getCostGradient(Vector& grad) { return getCostGradient(grad.getDataVector()); }

  static std::vector<ParamBound> getParameterBounds(int nt, int num_basis, bool linear_basis);

  void updateOptimalSolution(const double& cost_rel, const std::array<double,3>& sub_costs,
                             const Vector& param_vec);

  inline double uniformRand(double min = 0.0, double max = 1.0)
  {
    return min + (max - min)*uniform_rand(rand_engine);
  }
};

#endif
