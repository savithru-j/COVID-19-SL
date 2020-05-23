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


int main();

std::vector<Population> predictModel(const ModelParams& params, const Population& pop_init);

ModelParams getOptimalParameters(const ObservedPopulation& observed_pop);

void copyParam2Vector(const ModelParams& params, Vector& v);
void copyVector2Param(const Vector& v, ModelParams& params);

std::vector<ParamBound> getParameterBounds(int nt);

double uniformRand(double min = 0.0, double max = 1.0);

void randomizeParameterVector(const std::vector<ParamBound>& bounds, const double energy, Vector& param_vec);

double getOptCost(const ModelParams& params, const Vector& param_vec, const Population& pop_init,
                  const ObservedPopulation& observed_pop, const Matrix& regularization_matrix = Matrix(0,0),
                  const double scaling = 1.0);

Matrix getHaarMatrix(int m);

#endif
