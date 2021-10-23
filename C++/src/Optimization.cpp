//#include <random>

#include "Optimization.h"
#include "LinearAlgebra.h"

//double uniformRand(double min, double max)
//{
////  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
//  static std::mt19937 gen(2.0); //Standard mersenne_twister_engine seeded with rd()
//  static std::uniform_real_distribution<> uniform_rand(0, 1);
//  return min + (max - min)*uniform_rand(gen);
//}

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

std::string
getNLOPTResultDescription(nlopt::result resultcode)
{
  switch (resultcode)
  {
    case nlopt::result::FAILURE:
      return "Failed - Generic result code";
      break;

    case nlopt::result::INVALID_ARGS:
      return "Failed - Invalid arguments";
      break;

    case nlopt::result::OUT_OF_MEMORY:
      return "Failed - Out of memory";
      break;

    case nlopt::result::ROUNDOFF_LIMITED:
      return "Failed - Round-off limited";
      break;

    case nlopt::result::FORCED_STOP:
      return "Failed - Forcefully stopped";
      break;

    case nlopt::result::SUCCESS:
      return "Success - Generic result code";
      break;

    case nlopt::result::STOPVAL_REACHED:
      return "Success - Stop value reached";
      break;

    case nlopt::result::FTOL_REACHED:
      return "Success - Relative f-tolerance reached";
      break;

    case nlopt::result::XTOL_REACHED:
      return "Success - Relative x-tolerance reached";
      break;

    case nlopt::result::MAXEVAL_REACHED:
      return "Success - Maximum evaluation count reached";
      break;

    case nlopt::result::MAXTIME_REACHED:
      return "Success - Maximum time reached";
      break;

    default:
      break;
  }

  return "Unknown result code";
}
