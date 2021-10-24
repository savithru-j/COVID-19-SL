#ifndef MODELPARAMS_H
#define MODELPARAMS_H

#include "LinearAlgebra.h"

template<class T = double>
struct ModelParams
{
  ModelParams(int nt_hist_, int nt_pred_, const Vector<double>& external_quarantine_input = {},
              const Vector<double>& vaccination_data = {},
              double betaN_val = 0.3, double ce_val = 0.0, double c0_val = 0.0,
              double c1_val = 0.0, double c2_val = 1.0, double c3_val = 1.0)
    : nt_hist(nt_hist_), nt_pred(nt_pred_)
  {
    const int nt = nt_hist + nt_pred;
    quarantine_input.resize(nt, 0);
    daily_vaccinations.resize(nt, 0);
    betaN.resize(nt);
    ce.resize(nt);
    c0.resize(nt);
    c1.resize(nt);
    c2.resize(nt);
    c3.resize(nt);
    IFR.resize(nt);

    for (int i = 0; i < nt; ++i)
    {
      betaN[i] = betaN_val;
      ce[i]    = ce_val;
      c0[i]    = c0_val;
      c1[i]    = c1_val;
      c2[i]    = c2_val;
      c3[i]    = c3_val;
      IFR[i]   = 0.005; //default value of 0.5%

      if (i < (int) external_quarantine_input.size())
        quarantine_input[i] = external_quarantine_input[i];

      if (i < (int) vaccination_data.size())
        daily_vaccinations[i] = vaccination_data[i];
    }
  }

  template<class T2>
  ModelParams(const ModelParams<T2>& other) :
    nt_hist(other.nt_hist), nt_pred(other.nt_pred),
    quarantine_input(other.quarantine_input), daily_vaccinations(other.daily_vaccinations),
    T_incub0(other.T_incub0), T_incub1(other.T_incub1), T_asympt(other.T_asympt),
    T_mild(other.T_mild), T_severe(other.T_severe), T_icu(other.T_icu), f(other.f),
    frac_recover_I1(other.frac_recover_I1), frac_recover_I2(other.frac_recover_I2),
    S_Reff(other.S_Reff), vaccine_eff(other.vaccine_eff)
  {
    betaN.resize(other.betaN.size());
    std::copy(other.betaN.cbegin(), other.betaN.cend(), betaN.begin());

    ce.resize(other.ce.size());
    std::copy(other.ce.cbegin(), other.ce.cend(), ce.begin());

    c0.resize(other.c0.size());
    std::copy(other.c0.cbegin(), other.c0.cend(), c0.begin());

    c1.resize(other.c1.size());
    std::copy(other.c1.cbegin(), other.c1.cend(), c1.begin());

    c2.resize(other.c2.size());
    std::copy(other.c2.cbegin(), other.c2.cend(), c2.begin());

    c3.resize(other.c3.size());
    std::copy(other.c3.cbegin(), other.c3.cend(), c3.begin());

    IFR.resize(other.IFR.size());
    std::copy(other.IFR.cbegin(), other.IFR.cend(), IFR.begin());
  }

  int nt_hist, nt_pred;
  static constexpr double dt = 1.0/24.0;
  Vector<double> quarantine_input;
  Vector<double> daily_vaccinations;
  Vector<T> betaN;
  Vector<T> ce;
  Vector<T> c0;
  Vector<T> c1;
  Vector<T> c2;
  Vector<T> c3;
  double T_incub0        = 3.0;
  double T_incub1        = 2.0;
  double T_asympt        = 6.0;
  double T_mild          = 6.0;
  double T_severe        = 4.0;
  double T_icu           = 10.0;
  double f               = 0.3;   //exposed to asymptomatic probability
  double frac_recover_I1 = 0.80;  //fraction of cases that recover from mild-infected stage I1
  double frac_recover_I2 = 0.75;  //fraction of cases that recover from severe-infected stage I2
  Vector<T> IFR;                  //infection fatality ratio (for each day)
  double S_Reff          = 0.0;   //susceptible population at a given time - only used for R-effective calc
  T vaccine_eff          = 1.0;   //effectiveness of vaccines
};

template<class T>
void copyParam2FullVector(const ModelParams<T>& params, Vector<T>& v);

template<class T>
void copyFullVector2Param(const Vector<T>& v, ModelParams<T>& params);

#endif
