// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var block_size = 3;

function getParameterVector(params)
{
  return getParameterVector1(params);
}

function getParameterVector1(params)
{
  let Nt = params.T_hist - 1;
  let num_blocks = Math.ceil(Nt / block_size);

  let param_vec = [];

  for (let i = 0; i < num_blocks; ++i)
  {
    let sum = 0;
    let cnt = 0;
    for (let j = 0; j < block_size; ++j)
    {
      let ind = block_size*i + j;
      if (ind < Nt)
      {
        sum += params.b1N[ind];
        cnt++;
      }
    }
    param_vec.push(sum/cnt);
  }

  for (let i = 0; i < num_blocks; ++i)
  {
    let sum = 0;
    let cnt = 0;
    for (let j = 0; j < block_size; ++j)
    {
      let ind = block_size*i + j;
      if (ind < Nt)
      {
        sum += params.diag_frac[ind];
        cnt++;
      }
    }
    param_vec.push(sum/cnt);
  }

  // param_vec.push(params.E0_0);
  // param_vec.push(params.a0);
  // param_vec.push(params.a10);
  // param_vec.push(params.a11);
  // param_vec.push(params.g0);
  // param_vec.push(params.g1);
  // param_vec.push(params.p1);
  // param_vec.push(params.g2);
  // param_vec.push(params.p2);
  // param_vec.push(params.g3);
  // param_vec.push(params.mu);
  return param_vec;
}

function getParameterVector2(params)
{
  const T_end = params.T_hist - 1;
  const Nseg = 12;

  const n = (Nseg-1) + 2*Nseg;
  let param_vec = new Array(n).fill(0);

  let ind_vec = [0];
  for (let i = 0; i < Nseg; ++i)
    ind_vec.push(Math.round(T_end*(i+1)/Nseg));

  for (let i = 0; i < Nseg; ++i)
  {
    let jmin = ind_vec[i];
    let jmax = ind_vec[i+1];

    if (i < Nseg-1)
      param_vec[i] = jmax;

    let avg_b1N = 0, avg_c = 0;
    for (let j = jmin; j < jmax; ++j)
    {
      avg_b1N += params.b1N[j];
      avg_c += params.diag_frac[j];
    }
    avg_b1N /= (jmax - jmin);
    avg_c /= (jmax - jmin);

    param_vec[Nseg-1 + i] = avg_b1N;
    param_vec[Nseg-1 + Nseg + i] = avg_c;
  }

  return param_vec;
}

function updateParameterStructFromVector(params, param_vec, extrapolate_to_end = false)
{
  updateParameterStructFromVector1(params, param_vec);

  if (extrapolate_to_end)
  {
    let i_last = (params.T_hist - 1);
    for (let i = i_last; i < params.b1N.length; ++i)
      params.b1N[i] = params.b1N[i_last-1];

    for (let i = i_last; i < params.diag_frac.length; ++i)
      params.diag_frac[i] = params.diag_frac[i_last-1];
  }
}

function updateParameterStructFromVector1(params, param_vec)
{
  let Nt = params.T_hist - 1;
  let num_blocks = Math.ceil(Nt / block_size);

  for (let i = 0; i < num_blocks; ++i)
  {
    for (let j = 0; j < block_size; ++j)
    {
      let ind = block_size*i + j;
      if (ind < Nt)
      {
        params.b1N[ind] = param_vec[i];
        params.diag_frac[ind] = param_vec[num_blocks + i]
      }
    }
  }
  // let off = 2*num_blocks;
  // params.E0_0 = param_vec[off];
  // params.a0   = param_vec[off + 1];
  // params.a10  = param_vec[off + 2];
  // params.a11  = param_vec[off + 3];
  // params.g0   = param_vec[off + 4];
  // params.g1   = param_vec[off + 5];
  // params.p1   = param_vec[off + 6];
  // params.g2   = param_vec[off + 7];
  // params.p2   = param_vec[off + 8];
  // params.g3   = param_vec[off + 9];
  // params.mu   = param_vec[off + 10];
}

function updateParameterStructFromVector2(params, param_vec)
{
  const T_end = params.T_hist - 1;
  const Nseg = 12;

  let ind_vec = [0];
  for (let i = 0; i < Nseg-1; ++i)
    ind_vec.push(Math.round(param_vec[i]));
  ind_vec.push(T_end);

  for (let i = 0; i < Nseg; ++i)
  {
    let jmin = ind_vec[i];
    let jmax = ind_vec[i+1];

    for (let j = jmin; j < jmax; ++j)
    {
      params.b1N[j] = param_vec[Nseg-1 + i];
      params.diag_frac[j] = param_vec[Nseg-1 + Nseg + i];
    }
  }
}

function getParameterBounds(params)
{
  return getParameterBounds1(params);
}

function getParameterBounds1(params)
{
  let Nt = params.T_hist - 1;
  let num_blocks = Math.ceil(Nt / block_size);

  let bounds = [];

  for (let i = 0; i < num_blocks; ++i)
    bounds.push({min: 0.05, max: 1.0, step: 1e-4}); //b1N

  for (let i = 0; i < num_blocks; ++i)
    bounds.push({min: 0.05, max: 1.0, step: 1e-4}); //c

  // bounds.push({min: 1.0, max: 5.0, step: 1e-4}); //E0_0
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a0
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a10
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a11
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //g0
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //g1
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //p1
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //g2
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //p2
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //g3
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //mu

  return bounds;
}

function getLimitingStepSize(param_vec, dparam_vec, bounds)
{
  let eta = 1.0;

  for (let i = 0; i < param_vec.length; ++i)
  {
    let eta0 = eta;
    let new_val = param_vec[i] - dparam_vec[i];
    if (new_val < bounds[i].min)
      eta = Math.min(eta, (param_vec[i] - bounds[i].min)/dparam_vec[i]);
    else if (new_val > bounds[i].max)
      eta = Math.min(eta, (param_vec[i] - bounds[i].max)/dparam_vec[i]);

    if (eta < 0)
    {
      console.log("negative eta: " + eta + ", " + param_vec[i] + ", " + dparam_vec[i])
      break;
    }
  }
  return eta;
}

function randomizeParameterVector(param_vec, bounds)
{
  // for (let i = 0; i < param_vec.length; ++i)
  //   param_vec[i] = bounds[i].min + Math.random()*(bounds[i].max - bounds[i].min);

  for (let i = 0; i < param_vec.length; ++i)
  {
    let step = 2*Math.random() - 1.0;
    param_vec[i] += 0.1*step*(bounds[i].max - bounds[i].min);
    param_vec[i] = Math.min(Math.max(param_vec[i], bounds[i].min), bounds[i].max);
  }
}

function optimizeParameters()
{
  let T_pred_orig = sim_params.T_pred;
  let params = sim_params;
  params.T_pred = 0;

  //let params = initializeSimulationParameters(data_real.total.length, 0); //no prediction
  let param_vec = getParameterVector(params);
  let param_bounds = getParameterBounds(params);

  const n = param_vec.length;
  let dparam_vec = new Array(n).fill(0);
  let param_vec0 = new Array(n).fill(0);
  let param_vec_opt = new Array(n).fill(0);
  copyVector(param_vec, param_vec_opt);

  //Calculate initial cost
  let cost = getOptCost(params);
  let cost_init = cost;

  let cost_rel_opt = 1.0;
  const cost_reduction_tol = 1e-8;
  const min_eta = 1e-10;

  for (let pass = 0; pass < 50; ++pass)
  {
    cost = getOptCost(params);

    let cost_prev = cost;
    let stalled_iter = 0;

    let iter = 0;
    for (iter = 0; iter < 1000; ++iter)
    {
      let cost_grad = getOptCostGradient(params, param_bounds);

      //reset params struct, since Jacobian calc can change values slightly (by machine tolerances)
      updateParameterStructFromVector(params, param_vec);

      //Update parameter vector using gradient descent:
      //u(n+1) = u(n) - eta * (dR/du)^T R(u(n))
      for (let i = 0; i < n; ++i)
        param_vec0[i] = param_vec[i];

      let eta = 0.5*getLimitingStepSize(param_vec0, cost_grad, param_bounds);

      console.log("  Iter " + iter + ": " + cost.toExponential(4) + ", eta: " + eta.toExponential(4));

      if (eta < min_eta)
        break;

      while (eta >= min_eta)
      {
        //Update param_vec based on (scaled) update vector
        for (let i = 0; i < n; ++i)
          param_vec[i] = param_vec0[i] - eta*cost_grad[i];

        //This checks for the validity of parameters, and may modify values
        updateParameterStructFromVector(params, param_vec);

        //Evaluate new cost
        let cost_new = getOptCost(params);
        // console.log("  " + eta + ", " + getL2Norm(res));

        if (cost_new < cost)
        {
          cost = cost_new;
          break;
        }

        eta /= 2.0;
      } //linesearch

      //Check if residual has stalled
      if ((cost_prev - cost) < cost_prev*cost_reduction_tol)
        stalled_iter++;
      else {
        cost_prev = cost;
        stalled_iter = 0;
      }

      if (stalled_iter == 5) //Abort if residual has stalled for too many iterations.
        break;

    } //gradient descent

    let cost_rel = cost/cost_init;
    console.log("Pass " + pass + ", num_iter: " + iter + ", resnorm: " + cost.toExponential(4) + ", relative: " + cost_rel.toFixed(4));

    if (cost_rel < cost_rel_opt)
    {
      cost_rel_opt = cost_rel;
      copyVector(param_vec, param_vec_opt);
    }

    randomizeParameterVector(param_vec, param_bounds);
    updateParameterStructFromVector(params, param_vec);

  } //passes

  console.log("Optimal relative cost: " + cost_rel_opt.toFixed(4));

  //Copy final solution to global simulation parameters
  updateParameterStructFromVector(sim_params, param_vec_opt, true);
  // limitParameters(sim_params, true);
  sim_params.T_pred = T_pred_orig;

  // console.log("Final gamma: " + params.g0 + ", " + params.g1);

  updateParameters(true);
  return params;
}

function getOptCost(params)
{
  let sol_hist = predictModel(params);

  let wA = 1.0, wR = 1.0, wF = 1.0;

  let cost = 0.0;
  for (let i = 1; i < sol_hist.length; ++i)
  {
    //Error in number of active patients
    let num_active_pred = sol_hist[i][4] + sol_hist[i][6] + sol_hist[i][7] + sol_hist[i][8]; //I1d + I1q + I2 + I3
    let num_active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];

    // let cA = (num_active_true == 0) ? 1 : (1.0/num_active_true);
    // let cR = (data_real.categorized[i].y[1] == 0) ? 1 : (1.0/data_real.categorized[i].y[1]);
    // let cF = (data_real.categorized[i].y[2] == 0) ? 1 : (1.0/data_real.categorized[i].y[2]);

    let err_active = (num_active_pred - num_active_true);
    let err_recov = (sol_hist[i][9] - data_real.categorized[i].y[1]);
    let err_fatal = (sol_hist[i][11] - data_real.categorized[i].y[2]);
    let err_AR = (num_active_pred + sol_hist[i][9] - num_active_true - data_real.categorized[i].y[1]);

    // cost += wA*err_active*err_active + wR*err_recov*err_recov + wF*err_fatal*err_fatal;
    cost += wA*err_AR*err_AR + wF*err_fatal*err_fatal;
  }

  if (isNaN(cost))
  {
    console.log("getOptimizationCost: found NaN");
    return {};
  }
  return cost;
}

function getOptCostGradient(params, param_bounds)
{
  let param_vec = getParameterVector(params);
  const n = param_vec.length;

  let grad = new Array(n).fill(0);

  for (let j = 0; j < n; ++j)
  {
    let delta = param_bounds[j].step;

    //Compute finite difference
    param_vec[j] += delta;
    updateParameterStructFromVector(params, param_vec);
    let fp = getOptCost(params);

    param_vec[j] -= 2*delta;
    updateParameterStructFromVector(params, param_vec);
    let fm = getOptCost(params);

    param_vec[j] += delta;

    grad[j] = (fp - fm)/(2*delta);
    if (isNaN(grad[j]))
    {
      console.log("getOptCostGradient: found NaN");
      return {};
    }
  }
  return grad;
}

function getCurrentPredictionError()
{
  let res_sq = 0.0, res0_sq = 0.0;
  for (let i = 1; i < data_real.total.length; ++i)
  {
    let active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];
    let active_pred = 0;
    for (let j = 0; j < data_predicted.cat_diag.length; ++j)
      active_pred += data_predicted.cat_diag[j][i].y;

    let err_a = active_pred - active_true;
    let err_r = data_predicted.cat_diag[5][i].y - data_real.categorized[i].y[1]; //recovered
    let err_d = data_predicted.cat_diag[6][i].y - data_real.categorized[i].y[2]; //fatal

    res_sq += err_a*err_a + err_r*err_r + err_d*err_d;
    res0_sq += active_true*active_true + data_real.categorized[i].y[1]*data_real.categorized[i].y[1] + data_real.categorized[i].y[2]*data_real.categorized[i].y[2];
  }
  return Math.sqrt(res_sq/res0_sq);
}

function getL2Norm(vec)
{
  let norm = 0.0;
  for (let i = 0; i < vec.length; ++i)
    norm += vec[i]*vec[i];
  return Math.sqrt(norm);
}

function copyVector(vfrom, vto)
{
  for (let i = 0; i < vfrom.length; ++i)
    vto[i] = vfrom[i];
}
