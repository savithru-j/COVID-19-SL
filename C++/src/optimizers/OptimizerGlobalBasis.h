#ifndef OPTIMIZERGLOBALBASIS_H
#define OPTIMIZERGLOBALBASIS_H

#include <random>

#include "Optimization.h"
#include "models/Population.h"
#include "models/ObservedPopulation.h"
#include "linearalgebra/Vector.h"
#include "linearalgebra/Matrix.h"

struct OptimizerGlobalBasis
{
  static constexpr int NUM_RESULTS = 40;  //no. of optimal results to store (best to worst)

  OptimizerGlobalBasis(const ObservedPopulation& pop_observed_, const Population<double>& pop_init_,
                       const Vector<double>& quarantine_input, int num_basis_,
                       double wconf_ = 1, double wrecov_ = 1, double wfatal_ = 1,
                       int max_iter_per_pass_ = 1000, int max_passes_ = 1, int seed = 1);

  inline int nDim() const { return param_vec.size(); };

  inline void randomizeParameters()
  {
    for (int i = 0; i < 4; ++i)
    {
      param_vec[i*num_basis] = uniformRand(param_bounds[i*num_basis].min, param_bounds[i*num_basis].max);
      for (int j = 1; j < num_basis; ++j)
        param_vec[i*num_basis + j] = 0.0;
    }

    const int off = 4*num_basis;
    for (int i = off; i < param_vec.m(); ++i)
      param_vec[i] = uniformRand(param_bounds[i].min, param_bounds[i].max);

    copyVector2Param(param_vec, params);
  }

  void optimizeParametersNLOPT();

  const ObservedPopulation& pop_observed;
  const Population<double>& pop_init;
  const int nt_opt, num_basis;
  double weight_conf, weight_recov, weight_fatal;
  int max_iter_per_pass = 1000;
  int max_passes = 1;

  std::vector<ParamBound> param_bounds;
  Vector<double> param_vec; //current solution vector
  ModelParams<double> params;

  int f_eval_count = 0;
  int nlopt_iter = 0;

  std::array<Vector<double>,NUM_RESULTS> optimal_param_vec;
  std::array<double,NUM_RESULTS> cost_min;
  std::array<std::array<double,3>,NUM_RESULTS> sub_costs_min;

  std::mt19937 rand_engine; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> uniform_rand;

  Matrix<double> transform_matrix; //coeff_vec = A * raw_vec
  Matrix<double> inv_transform_matrix; //raw_vec = A * coeff_vec

  static double getCostNLOPT(const std::vector<double>& x, std::vector<double>& grad, void* data);

  static void getConstraintsNLOPT(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data);

  void copyParam2Vector(const ModelParams<double>& params, Vector<double>& v);
  void copyVector2Param(const Vector<double>& v, ModelParams<double>& params);

protected:

  std::pair<double, std::array<double, 3>> getCost();

  double getCostGradient(std::vector<double>& grad);
  double getCostGradient(Vector<double>& grad) { return getCostGradient(grad.getDataVector()); }

  void getConstraints(double* constraints, double* grad);
  void getConstraintsDCT(double* constraints, double* grad);

  double getPenalty();

  static std::vector<ParamBound> getParameterBounds(int nt, int num_basis);

//  static void evaluateLegendrePolynomial(const int nbasis, const int nt, const double* coeff, Vector& params);

  static void getBasisLegendre(const double& x, Vector<double>& phi);

  void updateOptimalSolution(const double& cost_rel, const std::array<double,3>& sub_costs,
                             const Vector<double>& param_vec);

  inline double uniformRand(double min = 0.0, double max = 1.0)
  {
    return min + (max - min)*uniform_rand(rand_engine);
  }
};

#endif
