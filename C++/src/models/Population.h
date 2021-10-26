#ifndef POPULATION_H
#define POPULATION_H

#include <array>
#include <fstream>

#include "ModelParams.h"
#include "utils/ErrorHandler.h"

template<class T = double>
class Population
{
public:

  using Array2 = std::array<T,2>;

  Population(double N_, double E0_0 = 0.0, double I1d_0 = 0.0, double Rd_0 = 0.0, double Dd_0 = 0.0)
    : N(N_), S(N - E0_0 - I1d_0 - Rd_0 - Dd_0), E0(E0_0), E1({0,0}),
      I0({0,0}), I1({0, I1d_0}), I2({0,0}), I3({0,0}), R({0, Rd_0}), D({0,Dd_0}) {}

  template<class T2>
  Population(const Population<T2>& other) :
    N(other.N), S(other.S), E0(other.E0)
  {
    std::copy(other.E1.cbegin(), other.E1.cend(), E1.begin());
    std::copy(other.I0.cbegin(), other.I0.cend(), I0.begin());
    std::copy(other.I1.cbegin(), other.I1.cend(), I1.begin());
    std::copy(other.I2.cbegin(), other.I2.cend(), I2.begin());
    std::copy(other.I3.cbegin(), other.I3.cend(), I3.begin());
    std::copy(other.R.cbegin(), other.R.cend(), R.begin());
    std::copy(other.D.cbegin(), other.D.cend(), D.begin());
    dS_exit_Reff = other.dS_exit_Reff;
  }

  void evolve(const ModelParams<T>& params, int t);
  void report(const ModelParams<T>& params, int t);
  void vaccinate(const ModelParams<T>& params, int t);

  inline T getNumReported() const { return E1[1] + I0[1] + I1[1] + I2[1] + I3[1] + R[1] + D[1]; }
  inline T getNumActiveReported() const { return E1[1] + I0[1] + I1[1] + I2[1] + I3[1] + R[1]; }
  inline T getNumRecoveredReported() const { return R[1]; }
  inline T getNumFatalReported() const { return D[1]; }

  inline T getNumUnreported() const { return E1[0] + I0[0] + I1[0] + I2[0] + I3[0] + R[0] + D[0]; }
  inline T getNumInfectedUnreported() const { return I0[0] + I1[0] + I2[0] + I3[0]; }
  inline T getNumRecoveredUnreported() const { return R[0]; }
  inline T getNumFatalUnreported() const { return D[0]; }

  double N = 0;
  T S = 0;
  T E0 = 0;
  Array2 E1;
  Array2 I0;
  Array2 I1;
  Array2 I2;
  Array2 I3;
  Array2 R;
  Array2 D;
  T dS_exit_Reff = 0;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Population<T>& pop);

#endif
