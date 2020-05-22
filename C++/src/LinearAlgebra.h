#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
#include "ErrorHandler.h"

class Vector
{
public:

  Vector(int m, double val = 0.0) : data_(m, val) {}

  Vector(const std::vector<double>& v) : data_(v) {}

  int m() const { return data_.size(); }

  const double& operator()(int i) const { return data_[i]; }
  double& operator()(int i) { return data_[i]; }

  const double& operator[](int i) const { return data_[i]; }
  double& operator[](int i) { return data_[i]; }

protected:
  std::vector<double> data_;
};

class Matrix
{
public:

  Matrix(int m, int n, double val = 0.0)
    : m_(m), n_(n), data_(m*n, val) {}

  int m() const { return m_; }
  int n() const { return n_; }

  const double& operator()(int i, int j) const { return data_[n_*i + j]; }
  double& operator()(int i, int j) { return data_[n_*i + j]; }

protected:
  int m_, n_;
  std::vector<double> data_;
};

Vector operator*(const Matrix& A, const Vector& x)
{
  if (A.n() != x.m())
    throwError("Dimension mismatch for matrix vector multiplication!");

  int m = A.m();
  int n = A.n();
  Vector y(m, 0.0);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      y(i) += A(i,j)*x(j);

  return y;
}

#endif
