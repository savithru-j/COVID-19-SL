#ifndef MODEL4LAYERPARAMS_H
#define MODEL4LAYERPARAMS_H

#include "linearalgebra/Vector.h"

template<class T = double>
struct ModelParams4Layer
{
  ModelParams4Layer(int nt_hist_, int nt_pred_, const Vector<double>& vaccination_data = {},
                    double beta_val = 0.3, double c_val = 0.0)
    : nt_hist(nt_hist_), nt_pred(nt_pred_)
  {
    const int nt = nt_hist + nt_pred;
    daily_vaccinations.resize(nt, 0);
    beta.resize(nt);
    c.resize(nt);
    IFR.resize(nt);

    for (int i = 0; i < nt; ++i)
    {
      beta[i] = beta_val;
      c[i]    = c_val;
      IFR[i]   = 0.005; //default value of 0.5%

      if (i < (int) vaccination_data.size())
        daily_vaccinations[i] = vaccination_data[i];
    }
  }

  template<class T2>
  ModelParams4Layer(const ModelParams4Layer<T2>& other) :
    nt_hist(other.nt_hist), nt_pred(other.nt_pred),
    daily_vaccinations(other.daily_vaccinations),
    T_incub(other.T_incub), T_recov(other.T_recov),
    S_Reff(other.S_Reff), beta_vac_scaling(other.beta_vac_scaling), 
    vaccine_alpha(other.vaccine_alpha), IFR_vac_scaling(other.IFR_vac_scaling)
  {
    beta.resize(other.beta.size());
    std::copy(other.beta.cbegin(), other.beta.cend(), beta.begin());

    c.resize(other.c.size());
    std::copy(other.c.cbegin(), other.c.cend(), c.begin());

    IFR.resize(other.IFR.size());
    std::copy(other.IFR.cbegin(), other.IFR.cend(), IFR.begin());
  }

  int nt_hist, nt_pred;
  static constexpr double dt = 1.0/24.0;
  Vector<double> daily_vaccinations;
  Vector<T> beta;
  Vector<T> c;
  double T_incub     = 5.0;   //time period from exposed state to infected state [days]
  double T_recov     = 6.0;   //time period from infected state to recovered state [days]
  Vector<T> IFR;              //infection fatality ratio (for each day)
  double S_Reff      = 0.0;   //susceptible population at a given time - only used for R-effective calc
  T beta_vac_scaling = 1.0;   //scaling to obtain beta for vaccinated infected people
  T vaccine_alpha    = 1.0;
  T IFR_vac_scaling  = 0.1;
  
};

template<class T>
void copyParam2FullVector(const ModelParams4Layer<T>& params, Vector<T>& v);

template<class T>
void copyFullVector2Param(const Vector<T>& v, ModelParams4Layer<T>& params);

#endif
