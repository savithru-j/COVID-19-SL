#include <iostream>
#include <iomanip>

#include "Population.h"

void
Population::evolve(const ModelParams& params, int t)
{
  //Probabilities
  const double prob_E0_E1 = 1;                      //non-infectious exposed to infectious exposed
  const double prob_E1_I0 = params.f;               //exposed to asymptomatic
  const double prob_E1_I1 = 1 - prob_E1_I0;         //exposed to mild
  const double prob_I0_R  = 1;
  const double prob_I1_R  = params.frac_recover_I1; //mild to recovered //0.8
  const double prob_I1_I2 = 1 - prob_I1_R;          //mild to severe  //0.2
  const double prob_I2_R  = params.frac_recover_I2; //severe to recovered  //0.75
  const double prob_I2_I3 = 1 - prob_I2_R;          //severe to critical //0.25
  const double prob_I3_D  = (prob_I1_I2 == 0 || prob_I2_I3 == 0) ?
                            0.0 : std::min(params.CFR/(prob_I1_I2*prob_I2_I3), 1.0); //critical to dead //0.40
  const double prob_I3_R  = 1 - prob_I3_D;          //critical to recovered //0.6

  //set rate parameters [1/day]
  const double a0  = (1/params.T_incub0) * prob_E0_E1;
  const double a10 = (1/params.T_incub1) * prob_E1_I0;
  const double a11 = (1/params.T_incub1) * prob_E1_I1;
  const double g0  = (1/params.T_asympt) * prob_I0_R;
  const double g1  = (1/params.T_mild)   * prob_I1_R;
  const double p1  = (1/params.T_mild)   * prob_I1_I2;
  const double g2  = (1/params.T_severe) * prob_I2_R;
  const double p2  = (1/params.T_severe) * prob_I2_I3;
  const double g3  = (1/params.T_icu)    * prob_I3_R;
  const double mu  = (1/params.T_icu)    * prob_I3_D;
  const double gd =  (1/params.T_discharge);

  double a1 = a10 + a11;
  double b = params.betaN[t] / N;

  //Only individuals in layer 0 (unreported) contribute to dS
  double dS = -b*(E1[0] + I0[0] + I1[0] + I2[0] + I3[0]) * S;
  double dE0 = -dS - a0*E0;
  Array2 dE1 = {a0*E0, 0};
  Array2 dI0 = {0, 0};
  Array2 dI1 = {0, params.quarantine_input[t]};
  Array2 dI2 = {0, 0};
  Array2 dI3 = {0, 0};
  Array2 dR  = {0, -gd*R[1]};
  Array2 dD  = {0, 0};
  double dRd = gd*R[1];

  for (int d = 0; d < 2; ++d) //loop over layers
  {
    dE1[d] += -a1*E1[d];
    dI0[d] += a10*E1[d] - g0*I0[d];
    dI1[d] += a11*E1[d] - (g1 + p1)*I1[d];
    dI2[d] +=  p1*I1[d] - (g2 + p2)*I2[d];
    dI3[d] +=  p2*I2[d] - (g3 + mu)*I3[d];
    dR[d]  +=  g0*I0[d] + g1*I1[d] + g2*I2[d] + g3*I3[d];
    dD[d]  +=  mu*I3[d];
  }

  //Update states
  S  +=  dS * params.dt;
  E0 += dE0 * params.dt;
  for (int d = 0; d < 2; ++d)
  {
    E1[d] += dE1[d] * params.dt;
    I0[d] += dI0[d] * params.dt;
    I1[d] += dI1[d] * params.dt;
    I2[d] += dI2[d] * params.dt;
    I3[d] += dI3[d] * params.dt;
    R[d]  +=  dR[d] * params.dt;
    D[d]  +=  dD[d] * params.dt;
  }
  Rd += dRd * params.dt;

  if (E0 < 0.5)
    E0 = 0.0;
}

void
Population::report(const ModelParams& params, int t)
{
  double delta = E1[0] * params.ce[t];
  E1[0] -= delta;
  E1[1] += delta;

  delta = I0[0] * params.c0[t];
  I0[0] -= delta;
  I0[1] += delta;

  delta = I1[0] * params.c1[t];
  I1[0] -= delta;
  I1[1] += delta;

  delta = I2[0] * params.c2[t];
  I2[0] -= delta;
  I2[1] += delta;

  delta = I3[0] * params.c3[t];
  I3[0] -= delta;
  I3[1] += delta;
}

std::ostream& operator<<(std::ostream& os, const Population& pop)
{
  os << std::scientific << std::setprecision(3);
  return os << "S : " << pop.S << std::endl
            << "E0: " << pop.E0 << std::endl
            << "E1: " << pop.E1[0] << ", " << pop.E1[1] << std::endl
            << "I0: " << pop.I0[0] << ", " << pop.I0[1] << std::endl
            << "I1: " << pop.I1[0] << ", " << pop.I1[1] << std::endl
            << "I2: " << pop.I2[0] << ", " << pop.I2[1] << std::endl
            << "I3: " << pop.I3[0] << ", " << pop.I3[1] << std::endl
            << "R : " << pop.R[0] << ", " << pop.R[1] << std::endl
            << "D : " << pop.D[0] << ", " << pop.D[1] << std::endl;
}

Vector
getQuarantineInputVector(const std::string& filepath)
{
  std::ifstream in(filepath);
  if (!in)
    return Vector(); //return empty vector

  std::vector<int> indices;
  std::vector<double> values;

  int ind, max_ind = -1;
  while (in >> ind)
  {
    double val;
    in >> val;
    indices.push_back(ind);
    values.push_back(val);
    max_ind = std::max(max_ind, ind);
  }
  in.close();

  Vector vec(max_ind+1, 0.0);
  for (std::size_t i = 0; i < indices.size(); ++i)
    vec[indices[i]] = values[i];

  return vec;
}
