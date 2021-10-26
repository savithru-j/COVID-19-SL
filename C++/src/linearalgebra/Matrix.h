#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#include "utils/ErrorHandler.h"

template<class T = double>
class Matrix
{
public:

  Matrix() = default;
  Matrix(int m, int n, double val = 0.0)
    : m_(m), n_(n), data_(m*n, val) {}

  inline int m() const { return m_; }
  inline int n() const { return n_; }

  inline const T& operator()(int i, int j) const { return data_[n_*i + j]; }
  inline T& operator()(int i, int j) { return data_[n_*i + j]; }

protected:
  int m_ = 0, n_ = 0;
  std::vector<T> data_;
};

template<class T>
inline Vector<T> operator*(const Matrix<T>& A, const Vector<T>& x)
{
  if (A.n() != x.m())
    throwError("Dimension mismatch for matrix vector multiplication!");

  int m = A.m();
  int n = A.n();
  Vector<T> y(m, 0.0);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      y(i) += A(i,j)*x(j);

  return y;
}

Matrix<double> getHaarMatrix(int m);
Matrix<double> getDCTMatrix(int N, int num_coeff);
Matrix<double> getInverseDCTMatrix(int N, int num_coeff);

#endif
