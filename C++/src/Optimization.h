#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <nlopt.hpp>

template<class T>
class Matrix;

struct ParamBound {

  ParamBound() = default;
  ParamBound(double min_, double max_, double step_) : min(min_), max(max_), step(step_) {}

  double min = 0.0, max = 1.0, step = 1e-4;
};

//double uniformRand(double min = 0.0, double max = 1.0);

Matrix<double> getHaarMatrix(int m);
Matrix<double> getDCTMatrix(int N, int num_coeff);
Matrix<double> getInverseDCTMatrix(int N, int num_coeff);

std::string getNLOPTResultDescription(nlopt::result resultcode);

#endif
