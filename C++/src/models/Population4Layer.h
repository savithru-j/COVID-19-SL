#ifndef POPULATION4LAYER_H
#define POPULATION4LAYER_H

#include <array>
#include <fstream>

#include "ModelParams4Layer.h"
#include "utils/ErrorHandler.h"

template<class T = double>
class Population4Layer
{
public:

  using Array4 = std::array<T,4>;

  Population4Layer(double N_, double E_0 = 0.0, double Ir_0 = 0.0, double Rr_0 = 0.0, double Dr_0 = 0.0)
    : N(N_), S(N - E_0 - Ir_0 - Rr_0 - Dr_0), Sv(0), E(E_0), Ev(0),
      I({0, Ir_0, 0, 0}), R({0, Rr_0, 0, 0}), D({0, Dr_0, 0, 0}) {}

  template<class T2>
  Population4Layer(const Population4Layer<T2>& other) :
    N(other.N), S(other.S), Sv(other.Sv), E(other.E), Ev(other.Ev)
  {
    std::copy(other.I.cbegin(), other.I.cend(), I.begin());
    std::copy(other.R.cbegin(), other.R.cend(), R.begin());
    std::copy(other.D.cbegin(), other.D.cend(), D.begin());
    dS_exit_Reff = other.dS_exit_Reff;
  }

  void evolve(const ModelParams4Layer<T>& params, int t);
  void report(const ModelParams4Layer<T>& params, int t);
  void vaccinate(const ModelParams4Layer<T>& params, int t);

  inline T getNumReported() const { return I[1] + I[3] + R[1] + R[3] + D[1] + D[3]; }
  inline T getNumActiveReported() const { return I[1] + I[3] + R[1] + R[3]; }
  inline T getNumRecoveredReported() const { return R[1] + R[3]; }
  inline T getNumFatalReported() const { return D[1] + D[3]; }

  inline T getNumUnreported() const { return I[0] + I[2] + R[0] + R[2] + D[0] + D[2]; }
  inline T getNumInfectedUnreported() const { return I[0] + I[2]; }
  inline T getNumRecoveredUnreported() const { return R[0] + R[2]; }
  inline T getNumFatalUnreported() const { return D[0] + D[2]; }

  double N = 0;
  T S = 0, Sv = 0; //susceptible (unvaccinated, vaccinated)
  T E = 0, Ev = 0; //exposed (unvaccinated, vaccinated)
  Array4 I; //infected (unvaccinated+unreported, unvaccinated+reported, vaccinated+unreported, vaccinated+reported)
  Array4 R; //reported (unvaccinated+unreported, unvaccinated+reported, vaccinated+unreported, vaccinated+reported)
  Array4 D; //deaths (unvaccinated+unreported, unvaccinated+reported, vaccinated+unreported, vaccinated+reported)
  T dS_exit_Reff = 0;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Population4Layer<T>& pop);

#endif
