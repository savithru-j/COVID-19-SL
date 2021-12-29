#include <cassert>
#include "Simulator.h"
#include "Population.h"
#include "Population2Layer.h"
#include "Population4Layer.h"
#include "ModelParams.h"
#include "ModelParams2Layer.h"
#include "ModelParams4Layer.h"
#include "linearalgebra/SurrealS.h"

template<class T>
T
calcEffectiveReproductionRatio(const ModelParams<T>& params_orig,
                               const std::vector<Population<T>>& pop_hist_orig, const int t)
{
  assert( t >= 0 && t < params_orig.nt_hist + params_orig.nt_pred);

  const int T_pred_Reff = 60;
  const double E0_init = 10000;
  const int nt_sub = 1.0/params_orig.dt;

  ModelParams<T> params(0, T_pred_Reff, {}, {}, params_orig.betaN[t], params_orig.ce[t], params_orig.c0[t],
                        params_orig.c1[t], params_orig.c2[t], params_orig.c3[t]);
  params.T_incub0        = params_orig.T_incub0;
  params.T_incub1        = params_orig.T_incub1;
  params.T_asympt        = params_orig.T_asympt;
  params.T_mild          = params_orig.T_mild;
  params.T_severe        = params_orig.T_severe;
  params.T_icu           = params_orig.T_icu;
  params.f               = params_orig.f;
  params.frac_recover_I1 = params_orig.frac_recover_I1;
  params.frac_recover_I2 = params_orig.frac_recover_I2;
  params.IFR             = params_orig.IFR;
  params.S_Reff          = pop_hist_orig[t].S;

  Population<T> pop(pop_hist_orig[t].N, E0_init);
  pop.S = 0; //susceptible population should be set to zero for R-eff simulation

  for (int t = 0; t < T_pred_Reff-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);
    pop.report(params, t); //report once daily
  }

  return pop.dS_exit_Reff / E0_init;
}

template<class T>
T
calcEffectiveReproductionRatio(const ModelParams2Layer<T>& params_orig,
                               const std::vector<Population2Layer<T>>& pop_hist_orig, const int t)
{
  assert(t >= 0 && t < params_orig.nt_hist + params_orig.nt_pred);

  const int T_pred_Reff = 90;
  const double E_init = 10000;
  const int nt_sub = 1.0/params_orig.dt;

  // ModelParams2Layer<T> params(0, T_pred_Reff, {}, {}, params_orig.beta[t]);
  // params.T_incub          = params_orig.T_incub;
  // params.T_recov          = params_orig.T_recov;
  // params.IFR              = params_orig.IFR;
  // params.beta_vac_scaling = params_orig.beta_vac_scaling;
  // params.vaccine_alpha    = params_orig.vaccine_alpha;
  // params.IFR_vac_scaling  = params_orig.IFR_vac_scaling;
  // params.S_Reff           = pop_hist_orig[t].S;
  // params.Sv_Reff          = pop_hist_orig[t].Sv;

  ModelParams2Layer<T> params(params_orig, t, T_pred_Reff);

  Population2Layer<T> pop(pop_hist_orig[t].N, E_init);
  pop.S = pop_hist_orig[t].S; //susceptible population should be set to zero for R-eff simulation
  pop.Sv = pop_hist_orig[t].Sv;

  for (int t = 0; t < T_pred_Reff-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t, true);
    pop.report(params, t); //report once daily
    pop.vaccinate(params, t); //remove vaccinated individuals daily
  }

  return pop.dS_exit_Reff / E_init;
}

template<class T>
T
calcEffectiveReproductionRatio(const ModelParams4Layer<T>& params_orig,
                               const std::vector<Population4Layer<T>>& pop_hist_orig, const int t)
{
  assert( t >= 0 && t < params_orig.nt_hist + params_orig.nt_pred);

  const int T_pred_Reff = 60;
  const double E_init = 10000;
  const int nt_sub = 1.0/params_orig.dt;

  ModelParams4Layer<T> params(0, T_pred_Reff, {}, params_orig.beta[t], params_orig.c[t]);
  params.T_incub          = params_orig.T_incub;
  params.T_recov          = params_orig.T_recov;
  params.IFR              = params_orig.IFR;
  params.beta_vac_scaling = params_orig.beta_vac_scaling;
  params.vaccine_alpha    = params_orig.vaccine_alpha;
  params.IFR_vac_scaling  = params_orig.IFR_vac_scaling;
  params.S_Reff           = pop_hist_orig[t].S;
  params.Sv_Reff          = pop_hist_orig[t].Sv;

  Population4Layer<T> pop(pop_hist_orig[t].N, E_init);
  pop.S = 0; //susceptible population should be set to zero for R-eff simulation
  pop.Sv = 0;

  for (int t = 0; t < T_pred_Reff-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);
    pop.report(params, t, true); //report once daily
  }

  return pop.dS_exit_Reff / E_init;
}

//Explicit instantiations
// template std::vector<Population<double>> predictModel(const ModelParams<double>&, const Population<double>&);
// template std::vector<Population<SurrealS<1,double>>> predictModel(const ModelParams<SurrealS<1,double>>&, const Population<double>&);

// template std::vector<Population4Layer<double>> predictModel(const ModelParams4Layer<double>&, const Population4Layer<double>&);
// template std::vector<Population4Layer<SurrealS<1,double>>> predictModel(const ModelParams4Layer<SurrealS<1,double>>&, const Population4Layer<double>&);

template double calcEffectiveReproductionRatio(const ModelParams<double>&, const std::vector<Population<double>>&, const int);
template double calcEffectiveReproductionRatio(const ModelParams2Layer<double>&, const std::vector<Population2Layer<double>>&, const int);
template double calcEffectiveReproductionRatio(const ModelParams4Layer<double>&, const std::vector<Population4Layer<double>>&, const int);
