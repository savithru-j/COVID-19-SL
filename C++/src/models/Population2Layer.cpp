#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "Population2Layer.h"
#include "linearalgebra/SurrealS.h"

template<class T>
void
Population2Layer<T>::evolve(const ModelParams2Layer<T>& params, int t, bool Reff_calc)
{
  //Probabilities
  const double prob_E_I = 1.0;            //exposed to infected
  const T      prob_I_D  = params.IFR[t]; //infected to dead (unvaccinated)
  const T      prob_I_R  = 1 - prob_I_D;  //infected to recovered

  const T prob_Iv_Dv = params.IFR_vac_scaling*params.IFR[t]; //infected to dead (vaccinated)
  const T prob_Iv_Rv = 1 - prob_Iv_Dv;

  //set rate parameters [1/day]
  const T theta  = (1/params.T_incub) * prob_E_I;
  const T gamma  = (1/params.T_recov) * prob_I_R;
  const T mu     = (1/params.T_recov) * prob_I_D;
  const T gamma_v = (1/params.T_recov) * prob_Iv_Rv;
  const T mu_v    = (1/params.T_recov) * prob_Iv_Dv;

  T beta_N = params.beta[t] / N;
  T beta_v_N = params.beta_vac_scaling * beta_N;
  T beta_t_N = params.beta_test_scaling * beta_N;

  //Only infected individuals in layers 0 and 2 (unreported) contribute to dS
  T beta_I_N_eff = (beta_N*I + beta_t_N*Ir + beta_v_N*Iv + beta_t_N*beta_v_N*Ivr);
  T dS      = -beta_I_N_eff * S;
  T dSv     = -params.vaccine_alpha*beta_I_N_eff * Sv;
  T dS_exit = -(dS + dSv);

  if (Reff_calc)
  {
    dS = 0;
    dSv = 0;
  }

  T dE   = -dS  - theta*E;
  T dEv  = -dSv - theta*Ev;
  T dI   = theta*E  - (gamma   + mu  )*I;
  T dIr  =          - (gamma   + mu  )*Ir;
  T dIv  = theta*Ev - (gamma_v + mu_v)*Iv;
  T dIvr =          - (gamma_v + mu_v)*Ivr;
  T dR   = gamma   * (I + Ir);
  T dRv  = gamma_v * (Iv + Ivr);
  T dD   = mu      * (I + Ir);
  T dDv  = mu_v    * (Iv + Ivr);

  // using namespace std;
  // if (isnan(dS) || isnan(dSv) || isnan(dE) || isnan(dEv))
  //   std::cout << "Found nan!" << std::endl;

  //Update states
  S   += dS   * params.dt;
  Sv  += dSv  * params.dt;
  E   += dE   * params.dt;
  Ev  += dEv  * params.dt;
  I   += dI   * params.dt;
  Ir  += dIr  * params.dt;
  Iv  += dIv  * params.dt;
  Ivr += dIvr * params.dt;
  R   += dR   * params.dt;
  Rv  += dRv  * params.dt;
  D   += dD   * params.dt;
  Dv  += dDv  * params.dt;
  dS_exit_Reff += dS_exit * params.dt;

  if (E < 0.1)
    E = 0.0;
  if (Ev < 0.1)
    Ev = 0.0;
}

template<class T>
void
Population2Layer<T>::report(const ModelParams2Layer<T>& params, int t)
{
  T num_reported_total = (t < (int) params.daily_tests.size()) ? params.daily_tests[t] : 0.0;
  T num_reportable = S + Sv + E + Ev + I + Iv;
  if (num_reportable == 0.0 || num_reported_total == 0.0)
    return;

  T delta = std::min((I/num_reportable) * num_reported_total, I);
  I  -= delta;
  Ir += delta;
  I_reported_cumulative += delta;

  delta = std::min((Iv/num_reportable) * num_reported_total, Iv);
  Iv  -= delta;
  Ivr += delta;
  I_reported_cumulative += delta;
}

template<class T>
void
Population2Layer<T>::vaccinate(const ModelParams2Layer<T>& params, int t)
{
  T num_active_unvac = S + E + I + Ir + R;
  if (num_active_unvac == 0.0 || params.daily_vaccinations[t] == 0.0)
    return;

  T num_vaccinated_S  = std::min((S / num_active_unvac) * params.daily_vaccinations[t], S);
  T num_vaccinated_E  = std::min((E / num_active_unvac) * params.daily_vaccinations[t], E);
  T num_vaccinated_I  = std::min((I / num_active_unvac) * params.daily_vaccinations[t], I);
  T num_vaccinated_Ir = std::min((Ir/ num_active_unvac) * params.daily_vaccinations[t], Ir);
  T num_vaccinated_R  = std::min((R / num_active_unvac) * params.daily_vaccinations[t], R);

  S  -= num_vaccinated_S;
  Sv += num_vaccinated_S;

  E  -= num_vaccinated_E;
  Ev += num_vaccinated_E;

  I  -= num_vaccinated_I;
  Iv += num_vaccinated_I;

  Ir  -= num_vaccinated_Ir;
  Ivr += num_vaccinated_Ir;

  R  -= num_vaccinated_R;
  Rv += num_vaccinated_R;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const Population2Layer<T>& pop)
{
  os << std::scientific << std::setprecision(3);
  return os << "S : " << pop.S << ", " << pop.Sv << std::endl
            << "E : " << pop.E << ", " << pop.Ev << std::endl
            << "I : " << pop.I << ", " << pop.Ir << ", " << pop.Iv << ", " << pop.Ivr << std::endl
            << "R : " << pop.R << ", " << pop.Rv << std::endl
            << "D : " << pop.D << ", " << pop.Dv << std::endl;
}

//Explicit instantiations
template class Population2Layer<double>;
template class Population2Layer<SurrealS<1,double>>;
