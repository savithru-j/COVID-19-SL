#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "Population4Layer.h"
#include "linearalgebra/SurrealS.h"

template<class T>
void
Population4Layer<T>::evolve(const ModelParams4Layer<T>& params, int t)
{
  //Probabilities
  const double prob_E_I = 1.0;            //exposed to infected
  const T      prob_I_D  = params.IFR[t]; //infected to dead (unvaccinated)
  const T      prob_I_R  = 1 - prob_I_D;  //infected to recovered

  const T prob_Iv_Dv = params.IFR_vac_scaling*params.IFR[t]; //infected to dead (vaccinated)
  const T prob_Iv_Rv = 1 - prob_Iv_Dv;

  //set rate parameters [1/day]
  const double theta  = (1/params.T_incub) * prob_E_I;
  const T gamma  = (1/params.T_recov) * prob_I_R;
  const T mu     = (1/params.T_recov) * prob_I_D;
  const T gamma_v = (1/params.T_recov) * prob_Iv_Rv;
  const T mu_v    = (1/params.T_recov) * prob_Iv_Dv;

  T beta_N = params.beta[t] / N;
  T beta_v_N = params.beta_vac_scaling * beta_N;

  //Only infected individuals in layers 0 and 2 (unreported) contribute to dS
  T dS      = -(beta_N*I[0] + beta_v_N*I[2]) * S;
  T dSv     = -params.vaccine_alpha*(beta_N*I[0] + beta_v_N*I[2]) * Sv;
  T dS_exit =  (beta_N*I[0] + beta_v_N*I[2]) * params.S_Reff;

  T dE      = -dS  - theta*E;
  T dEv     = -dSv - theta*Ev;
  Array4 dI = {theta*E  - (gamma   + mu  )*I[0], 
                        - (gamma   + mu  )*I[1],
               theta*Ev - (gamma_v + mu_v)*I[2], 
                        - (gamma_v + mu_v)*I[3]};
  Array4 dR = {gamma  *I[0], gamma  *I[1], 
               gamma_v*I[2], gamma_v*I[3]};
  Array4 dD = {mu  *I[0], mu  *I[1], 
               mu_v*I[2], mu_v*I[3]};

  //Update states
  S  += dS  * params.dt;
  Sv += dSv * params.dt;
  dS_exit_Reff += dS_exit * params.dt;
  E  += dE  * params.dt;
  Ev += dEv * params.dt;
  for (int d = 0; d < 4; ++d)
  {
    I[d] += dI[d] * params.dt;
    R[d] += dR[d] * params.dt;
    D[d] += dD[d] * params.dt;
  }

  if (E < 0.1)
    E = 0.0;
  if (Ev < 0.1)
    Ev = 0.0;
}

template<class T>
void
Population4Layer<T>::report(const ModelParams4Layer<T>& params, int t)
{
  T num_reported_total = params.c[t] * N;
  T num_reportable = S + E + I[0] + I[2];

  T delta = std::min((I[0]/num_reportable) * num_reported_total, I[0]);
  I[0] -= delta;
  I[1] += delta;

  delta = std::min((I[2]/num_reportable) * num_reported_total, I[2]);
  I[2] -= delta;
  I[3] += delta;
}

template<class T>
void
Population4Layer<T>::vaccinate(const ModelParams4Layer<T>& params, int t)
{
  T num_active_unvac = S + E + I[0] + I[1] + R[0] + R[1];
  T num_vaccinated_S  = std::min((S   /num_active_unvac) * params.daily_vaccinations[t], S);
  T num_vaccinated_E  = std::min((E   /num_active_unvac) * params.daily_vaccinations[t], E);
  T num_vaccinated_Iu = std::min((I[0]/num_active_unvac) * params.daily_vaccinations[t], I[0]);
  T num_vaccinated_Ir = std::min((I[1]/num_active_unvac) * params.daily_vaccinations[t], I[1]);
  T num_vaccinated_Ru = std::min((R[0]/num_active_unvac) * params.daily_vaccinations[t], R[0]);
  T num_vaccinated_Rr = std::min((R[1]/num_active_unvac) * params.daily_vaccinations[t], R[1]);

  S  -= num_vaccinated_S;
  Sv += num_vaccinated_S;

  E  -= num_vaccinated_E;
  Ev += num_vaccinated_E;

  I[0] -= num_vaccinated_Iu;
  I[2] += num_vaccinated_Iu;

  I[1] -= num_vaccinated_Ir;
  I[3] += num_vaccinated_Ir;

  R[0] -= num_vaccinated_Ru;
  R[2] += num_vaccinated_Ru;
  
  R[1] -= num_vaccinated_Rr;
  R[3] += num_vaccinated_Rr;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const Population4Layer<T>& pop)
{
  os << std::scientific << std::setprecision(3);
  return os << "S : " << pop.S << std::endl
            << "Sv: " << pop.Sv << std::endl
            << "E : " << pop.E << std::endl
            << "Ev: " << pop.Ev << std::endl
            << "I : " << pop.I[0] << ", " << pop.I[1] << ", " << pop.I[2]<< ", " << pop.I[3] << std::endl
            << "R : " << pop.R[0] << ", " << pop.R[1] << ", " << pop.R[2]<< ", " << pop.R[3] << std::endl
            << "D : " << pop.D[0] << ", " << pop.D[1] << ", " << pop.D[2]<< ", " << pop.D[3] << std::endl;
}

//Explicit instantiations
template class Population4Layer<double>;
template class Population4Layer<SurrealS<1,double>>;
