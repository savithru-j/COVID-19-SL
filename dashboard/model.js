// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

class Population
{
  constructor(N)
  {
    this.N = N;
    this.S = N;
    this.E0 = 0;
    this.E1 = [0, 0];
    this.I0 = [0, 0];
    this.I1 = [0, 0];
    this.I2 = [0, 0];
    this.I3 = [0, 0];
    this.R = [0, 0];
    this.D = [0, 0];
  }

  clone()
  {
    let p = new Population(this.N);
    p.S = this.S;
    p.E0 = this.E0;
    p.E1 = this.E1.slice();
    p.I0 = this.I0.slice();
    p.I1 = this.I1.slice();
    p.I2 = this.I2.slice();
    p.I3 = this.I3.slice();
    p.R = this.R.slice();
    p.D = this.D.slice();
    return p;
  }

  getNumExposed() { return Math.round(this.E0) + Math.round(this.E1[0]) + Math.round(this.E1[1]); }
  getNumAsymptomatic() { return Math.round(this.I0[0]) + Math.round(this.I0[1]); }
  getNumMild() { return Math.round(this.I1[0]) + Math.round(this.I1[1]); }
  getNumSevere() { return Math.round(this.I2[0]) + Math.round(this.I2[1]); }
  getNumCritical() { return Math.round(this.I3[0]) + Math.round(this.I3[1]); }
  getNumRecovered() { return Math.round(this.R[0]) + Math.round(this.R[1]); }
  getNumFatal() { return Math.round(this.D[0]) + Math.round(this.D[1]); }

  getNumDiagnosed() {
    return Math.round(this.E1[1]) + Math.round(this.I0[1]) + Math.round(this.I1[1])
         + Math.round(this.I2[1]) + Math.round(this.I3[1]) + Math.round(this.R[1]) + Math.round(this.D[1]);
  }

  evolve(params, t)
  {
    const b1 = params.b1N[t] / this.N;
    const b2 = params.b2N[t] / this.N;
    const b3 = params.b3N[t] / this.N;

    const a0 = params.a0;
    const a10 = params.a10;
    const a11 = params.a11;
    const a1 = a10 + a11;
    const g0 = params.g0;
    const g1 = params.g1;
    const p1 = params.p1;
    const g2 = params.g2;
    const p2 = params.p2;
    const g3 = params.g3;
    const mu = params.mu;

    //Only individuals in layer 0 (unreported) contribute to dS
    let dS = -(b1*(this.E1[0] + this.I0[0] + this.I1[0]) + b2*this.I2[0] + b3*this.I3[0]) * this.S;
    let dE0 = -dS - a0*this.E0;
    let dE1 = [a0*this.E0, 0];
    let dI0 = [0, 0];
    let dI1 = [0, params.quarantine_input[t]];
    let dI2 = [0, 0];
    let dI3 = [0, 0];
    let dR = [0, 0];
    let dD = [0, 0];

    for (let d = 0; d < 2; ++d) //loop over layers
    {
      dE1[d] += -a1*this.E1[d];
      dI0[d] += a10*this.E1[d] - g0*this.I0[d];
      dI1[d] += a11*this.E1[d] - (g1 + p1)*this.I1[d];
      dI2[d] +=  p1*this.I1[d] - (g2 + p2)*this.I2[d];
      dI3[d] +=  p2*this.I2[d] - (g3 + mu)*this.I3[d];
      dR[d]  +=  g0*this.I0[d] + g1*this.I1[d] + g2*this.I2[d] + g3*this.I3[d];
      dD[d]  +=  mu*this.I3[d];
    }

    //Update states
    const dt = params.dt;
    this.S  += dS * dt;
    this.E0 += dE0 * dt;
    for (let d = 0; d < 2; ++d)
    {
      this.E1[d] += dE1[d] * dt;
      this.I0[d] += dI0[d] * dt;
      this.I1[d] += dI1[d] * dt;
      this.I2[d] += dI2[d] * dt;
      this.I3[d] += dI3[d] * dt;
      this.R[d]  +=  dR[d] * dt;
      this.D[d]  +=  dD[d] * dt;
    }

    if (this.E0 < 0.5)
      this.E0 = 0.0;
  }

  report(params, t)
  {
    let delta = this.E1[0] * params.ce[t];
    this.E1[0] -= delta;
    this.E1[1] += delta;

    delta = this.I0[0] * params.c0[t];
    this.I0[0] -= delta;
    this.I0[1] += delta;

    delta = this.I1[0] * params.c1[t];
    this.I1[0] -= delta;
    this.I1[1] += delta;

    delta = this.I2[0] * params.c2[t];
    this.I2[0] -= delta;
    this.I2[1] += delta;

    delta = this.I3[0] * params.c3[t];
    this.I3[0] -= delta;
    this.I3[1] += delta;
  }
}
