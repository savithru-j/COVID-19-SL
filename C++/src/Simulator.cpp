#include <cassert>
#include "Simulator.h"
#include "Population.h"
#include "ModelParams.h"

std::vector<Population>
predictModel(const ModelParams& params, const Population& pop_init)
{
  std::vector<Population> population_hist;

  const int nt = params.nt_hist + params.nt_pred;
  const int nt_sub = 1.0/params.dt;

  Population pop = pop_init;

  population_hist.reserve(nt);
  population_hist.push_back(pop);

  for (int t = 0; t < nt-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);

    pop.report(params, t); //report once daily
    pop.vaccinate(params, t); //remove vaccinated individuals daily

    population_hist.push_back(pop); //save solution daily
  }

  return population_hist;
}

double
calcEffectiveReproductionRatio(const ModelParams& params_orig,
                               const std::vector<Population>& pop_hist_orig, const int t)
{
  assert( t >= 0 && t < params_orig.nt_hist + params_orig.nt_pred);

  const int T_pred_Reff = 60;
  const double E0_init = 10000;
  const int nt_sub = 1.0/params_orig.dt;

  ModelParams params(0, T_pred_Reff, {}, params_orig.betaN[t], params_orig.ce[t], params_orig.c0[t],
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

  Population pop(pop_hist_orig[t].N, E0_init);
  pop.S = 0; //susceptible population should be set to zero for R-eff simulation

  for (int t = 0; t < T_pred_Reff-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);
    pop.report(params, t); //report once daily
  }

  return pop.dS_exit_Reff / E0_init;
}
