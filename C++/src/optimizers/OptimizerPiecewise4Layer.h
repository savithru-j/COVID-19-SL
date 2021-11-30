#ifndef OPTIMIZERPIECEWISE4LAYER_H
#define OPTIMIZERPIECEWISE4LAYER_H

#include <random>

#include "Optimization.h"
#include "models/Population4Layer.h"
#include "models/ObservedPopulation.h"
#include "linearalgebra/Vector.h"

struct OptimizerPiecewise4Layer
{
  // static constexpr bool OPTIMIZE_C0 = true;
  // static constexpr bool OPTIMIZE_C1 = true;
  // static constexpr bool OPTIMIZE_C2 = false;
  static constexpr std::array<int,3> IFR_SEG_STARTS = {0, 151, 273};

  OptimizerPiecewise4Layer(const ObservedPopulation& pop_observed_, const Population4Layer<double>& pop_init_,
                           const Vector<double>& vaccination_data,
                           int interval_size_, bool linear_basis_ = false,
                           double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1,
                           int max_iter_per_pass_ = 1000, int max_passes_ = 1, int seed = 1);

  inline int nDim() const { return param_bounds.size(); };

  inline void randomizeParameters(std::vector<double>& param_vec)
  {
    if (param_vec.size() != param_bounds.size())
      throwError("randomizeParameters - inconsistent dimensions!");

    const int num_nodes = (int)(nt_opt/interval_size) + 1*linear_basis;
    constexpr int num_dynamic_params = 2; //beta, c
    for (int i = 0; i < num_dynamic_params; ++i)
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
    int off = num_dynamic_params*num_nodes;
    param_vec[off] = uniformRand(param_bounds[off].min, param_bounds[off].max);
    for (std::size_t i = 1; i < IFR_SEG_STARTS.size(); ++i)
      param_vec[off + i] = param_vec[off];

    off += IFR_SEG_STARTS.size();
    for (std::size_t i = off; i < param_vec.size(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);
  }

  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population4Layer<double>& pop_init;
  const int nt_opt, interval_size;
  const bool linear_basis = false;
  double weight_conf, weight_recov, weight_fatal;
  int max_iter_per_pass = 1000;
  int max_passes = 1;

  Vector<ParamBound> param_bounds;
  ModelParams4Layer<double> params;

  int f_eval_count = 0;
  int nlopt_iter = 0;

  Vector<Vector<double>> optimal_param_vec;
  Vector<double> cost_min;
  Vector<std::array<double,3>> sub_costs_min;

  std::mt19937 rand_engine; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> uniform_rand;

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  template<class T>
  void copyParam2Vector(const ModelParams4Layer<T>& params, std::vector<T>& v);

  template<class T>
  inline void copyParam2Vector(const ModelParams4Layer<T>& params, Vector<T>& v) {
    copyParam2Vector(params, v.getDataVector());
  }

  template<class T>
  void copyVector2Param(const std::vector<T>& v, ModelParams4Layer<T>& params);

  template<class T>
  inline void copyVector2Param(const Vector<T>& v, ModelParams4Layer<T>& params) {
    copyVector2Param(v.getDataVector(), params);
  }

protected:

  template<class T>
  std::pair<T, std::array<T, 3>> getCost(const ModelParams4Layer<T>& params_tmp);

  //Evaluates the cost gradient using finite differences
  double getCostGradientFD(const ModelParams4Layer<double>& params, std::vector<double>& grad);

  //Evaluates the cost gradient using automatic differentiation
  double getCostGradientAD(const ModelParams4Layer<double>& params, std::vector<double>& grad);

  static Vector<ParamBound> getParameterBounds(int nt, int num_basis, bool linear_basis);

  void updateOptimalSolution(const double& cost_rel, const std::array<double,3>& sub_costs,
                             const Vector<double>& param_vec);

  inline double uniformRand(double min = 0.0, double max = 1.0)
  {
    return min + (max - min)*uniform_rand(rand_engine);
  }
};

#endif
