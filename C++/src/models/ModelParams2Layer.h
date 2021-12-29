#ifndef MODEL2LAYERPARAMS_H
#define MODEL2LAYERPARAMS_H

#include "linearalgebra/Vector.h"

template<class T = double>
struct ModelParams2Layer
{
  ModelParams2Layer(int nt_hist_, int nt_pred_, 
                    const Vector<double>& testing_data = {}, const Vector<double>& vaccination_data = {},
                    double beta_val = 0.3)
    : nt_hist(nt_hist_), nt_pred(nt_pred_)
  {
    const int nt = nt_hist + nt_pred;
    daily_tests.resize(nt, 0);
    daily_vaccinations.resize(nt, 0);
    beta.resize(nt);
    IFR.resize(nt);

    for (int i = 0; i < nt; ++i)
    {
      beta[i] = beta_val;
      IFR[i]  = 0.005; //default value of 0.5%

      if (i > 0 && i < (int) testing_data.size())
        daily_tests[i] = testing_data[i] - testing_data[i-1]; //Compute diff to get daily values

      if (i > 0 && i < (int) vaccination_data.size())
        daily_vaccinations[i] = vaccination_data[i] - vaccination_data[i-1]; //Compute diff to get daily values
    }
  }

  template<class T2>
  ModelParams2Layer(const ModelParams2Layer<T2>& other) :
    nt_hist(other.nt_hist), nt_pred(other.nt_pred),
    daily_tests(other.daily_tests),
    daily_vaccinations(other.daily_vaccinations),
    T_incub(other.T_incub), T_recov(other.T_recov),
    // S_Reff(other.S_Reff), Sv_Reff(other.Sv_Reff), 
    beta_test_scaling(other.beta_test_scaling),
    beta_vac_scaling(other.beta_vac_scaling),
    vaccine_alpha(other.vaccine_alpha), IFR_vac_scaling(other.IFR_vac_scaling)
  {
    beta.resize(other.beta.size());
    std::copy(other.beta.cbegin(), other.beta.cend(), beta.begin());

    IFR.resize(other.IFR.size());
    std::copy(other.IFR.cbegin(), other.IFR.cend(), IFR.begin());
  }

  template<class T2>
  ModelParams2Layer(const ModelParams2Layer<T2>& params_orig, const int t_cur, const int nt_pred_) :
    nt_hist(0), nt_pred(nt_pred_)
  {
    daily_tests.resize(nt_pred, 0);
    daily_vaccinations.resize(nt_pred, 0);
    beta.resize(nt_pred);
    IFR.resize(nt_pred);

    for (int i = 0; i < nt_pred; ++i)
    {
      beta[i] = params_orig.beta[t_cur];
      IFR[i]  = params_orig.IFR[t_cur];

      if (t_cur < (int) params_orig.daily_tests.size())
        daily_tests[i] = params_orig.daily_tests[t_cur];

      if (t_cur < (int) params_orig.daily_vaccinations.size())
        daily_vaccinations[i] = params_orig.daily_vaccinations[t_cur];
    }

    T_incub = params_orig.T_incub;
    T_recov = params_orig.T_recov;
    beta_test_scaling = params_orig.beta_test_scaling;
    beta_vac_scaling = params_orig.beta_vac_scaling;
    vaccine_alpha = params_orig.vaccine_alpha;
    IFR_vac_scaling = params_orig.IFR_vac_scaling;
  }

  int nt_hist, nt_pred;
  static constexpr double dt = 1.0/24.0;
  Vector<double> daily_tests;
  Vector<double> daily_vaccinations;
  Vector<T> beta;
  T T_incub          = 5.0;   //time period from exposed state to infected state [days]
  T T_recov          = 8.0;   //time period from infected state to recovered state [days]
  Vector<T> IFR;              //infection fatality ratio (for each day)
  // double S_Reff      = 0.0;   //susceptible population at a given time - only used for R-effective calc
  // double Sv_Reff     = 0.0;   //vaccinated susceptible population at a given time - only used for R-effective calc
  T beta_test_scaling = 0.5;  //scaling to obtain beta for reported infected people in isolation
  T beta_vac_scaling = 0.5;   //scaling to obtain beta for vaccinated infected people
  T vaccine_alpha    = 0.5;
  T IFR_vac_scaling  = 0.1;
};

template<class T>
void copyParam2FullVector(const ModelParams2Layer<T>& params, Vector<T>& v);

template<class T>
void copyFullVector2Param(const Vector<T>& v, ModelParams2Layer<T>& params);

#endif
