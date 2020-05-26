#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
#include "ErrorHandler.h"

class Vector
{
public:

  Vector() = default;
  Vector(int m, double val = 0.0) : data_(m, val) {}

  Vector(const std::vector<double>& v) : data_(v) {}

  int m() const { return data_.size(); }
  int size() const { return data_.size(); }

  const double& operator()(int i) const { return data_[i]; }
  double& operator()(int i) { return data_[i]; }

  const double& operator[](int i) const { return data_[i]; }
  double& operator[](int i) { return data_[i]; }

  const double& back() const { return data_.back(); }

  void resize(int m) { data_.resize(m); }

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

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
  os << "[";
  for (int i = 0; i < v.m()-1; ++i)
    os << v[i] << ", ";
  os << v.back() << "]";
  return os;
}

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& v)
{
  for (std::size_t i = 0; i < N-1; ++i)
    os << v[i] << ", ";
  os << v.back();
  return os;
}

#endif
