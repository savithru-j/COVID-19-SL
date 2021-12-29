#ifndef POPULATION2LAYER_H
#define POPULATION2LAYER_H

#include <array>
#include <fstream>

#include "ModelParams2Layer.h"
#include "utils/ErrorHandler.h"

template<class T = double>
class Population2Layer
{
public:

  Population2Layer(double N_, double E_0 = 0.0, double Ir_0 = 0.0, double Rr_0 = 0.0, double Dr_0 = 0.0)
    : N(N_), S(N - E_0 - Ir_0 - Rr_0 - Dr_0), Sv(0), 
      E(E_0), Ev(0),
      I(0), Ir(Ir_0), Iv(0), Ivr(0),
      R(Rr_0), Rv(0),
      D(Dr_0), Dv(0),
      I_reported_cumulative(0) {}

  template<class T2>
  Population2Layer(const Population2Layer<T2>& pop) :
    N(pop.N), S(pop.S), Sv(pop.Sv), E(pop.E), Ev(pop.Ev),
    I(pop.I), Ir(pop.Ir), Iv(pop.Iv), Ivr(pop.Ivr),
    R(pop.R), Rv(pop.Rv), D(pop.D), Dv(pop.Dv),
    dS_exit_Reff(pop.dS_exit_Reff), I_reported_cumulative(pop.I_reported_cumulative) {}

  void evolve(const ModelParams2Layer<T>& params, int t, bool Reff_calc = false);
  void report(const ModelParams2Layer<T>& params, int t);
  void vaccinate(const ModelParams2Layer<T>& params, int t);

  inline T getNumReported() const { return I_reported_cumulative; }
  inline T getNumInfectedReported() const { return Ir + Ivr; }
  // inline T getNumRecoveredReported() const { return R[1] + R[3]; }
  inline T getNumFatalReported() const { return D + Dv; }

  inline T getNumUnreported() const { return I + Ir + Iv + Ivr + R + Rv + D + Dv - I_reported_cumulative; }
  inline T getNumInfectedUnreported() const { return I + Iv; }
  // inline T getNumRecoveredUnreported() const { return R[0] + R[2]; }
  inline T getNumFatalUnreported() const { return 0; }

  inline T getNumVaccinated() const { return Sv + Ev + Iv + Ivr + Rv + Dv; }

  double N = 0;
  T S = 0, Sv = 0; //susceptible (unvaccinated, vaccinated)
  T E = 0, Ev = 0; //exposed (unvaccinated, vaccinated)
  T I = 0, Ir = 0, Iv = 0, Ivr = 0; //infected (unvaccinated+unreported, unvaccinated+reported, vaccinated+unreported, vaccinated+reported)
  T R = 0, Rv = 0; //reported (unvaccinated, vaccinated)
  T D = 0, Dv = 0; //deaths (unvaccinated, vaccinated)
  T dS_exit_Reff = 0;
  T I_reported_cumulative = 0; //unvaccinated + vaccinated
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Population2Layer<T>& pop);

#endif
