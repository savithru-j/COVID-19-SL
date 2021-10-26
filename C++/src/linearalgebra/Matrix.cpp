#include <cmath>
#include "Matrix.h"

Matrix<double> getHaarMatrix(int m)
{
  const int num_levels = std::ceil(std::log2(m));
  const int n = std::pow(2, num_levels);
//  if (std::pow(2,num_levels) != n)
//    throw("N is not a power of 2!");

  Matrix<double> A(n, n, 0.0);
  const double inv_sqrtn = 1.0 / std::sqrt(n);

  for (int j = 0; j < n; ++j)
    A(0, j) = inv_sqrtn;

  for (int k = 1; k < n; ++k)
  {
    int p = std::floor(std::log(k) / std::log(2.0));
    double k1 = std::pow(2.0, p);
    double k2 = k1 * 2.0;
    double q = k - k1;
    double t1 = n/k1;
    double t2 = n/k2;
    double tmp = std::pow(2.0, p/2.0) * inv_sqrtn;
    for (int i = 0; i < t2; ++i)
    {
      A(k, q*t1 + i) = tmp;
      A(k, + q*t1 + i + t2) = -tmp;
    }
  }

  return A;
}

Matrix<double> getDCTMatrix(int N, int num_coeff)
{
  if (num_coeff == 0)
    throw("getDCTMatrix: the number of coefficients cannot be zero!");

  const double scale_row0 = 1.0 / std::sqrt(N);
  const double scale_rowk = std::sqrt(2.0) / std::sqrt(N);

  Matrix<double> A(N, N, 0.0);

  //A(0,:)
  for (int n = 0; n < N; ++n)
    A(0,n) = scale_row0; //row k = 0

  //A(1:N,:)
  for (int k = 1; k < num_coeff; ++k)
    for (int n = 0; n < N; ++n)
      A(k,n) = std::cos(M_PI*(n + 0.5)*k / N) * scale_rowk;

  return A;
}

Matrix<double> getInverseDCTMatrix(int N, int num_coeff)
{
  if (num_coeff == 0)
    throw("getInverseDCTMatrix: the number of coefficients cannot be zero!");

  const double scale_col0 = 1.0 / std::sqrt(N);
  const double scale_colk = std::sqrt(2.0) / std::sqrt(N);

  Matrix<double> A(N, num_coeff, 0.0);

  for (int n = 0; n < N; ++n)
  {
    A(n,0) = scale_col0;

    for (int k = 1; k < num_coeff; ++k)
      A(n,k) = std::cos(M_PI*(n + 0.5)*k / N) * scale_colk;
  }
  return A;
}
