// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var block_size = 1;
var tail_block_size = 2;
const opt_array_names = ["b1N","c0","c1"];

function getParameterVector(params)
{
  return getParameterVector1(params);
}

function getParameterVector1(params)
{
  let Nt = params.T_hist - tail_block_size;
  let num_blocks = Math.ceil(Nt / block_size);

  let param_vec = [];

  for (let name of opt_array_names)
  {
    for (let i = 0; i < num_blocks; ++i)
    {
      let sum = 0;
      let cnt = 0;
      for (let j = 0; j < block_size; ++j)
      {
        let ind = block_size*i + j;
        if (ind < Nt)
        {
          sum += params[name][ind];
          cnt++;
        }
      }
      param_vec.push(sum/cnt);
    }
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

function updateParameterStructFromVector(params, param_vec, extrapolate_to_end = false)
{
  updateParameterStructFromVector1(params, param_vec);

  if (extrapolate_to_end)
  {
    let i_last = (params.T_hist - tail_block_size);
    for (let name of opt_array_names)
    {
      let val = params[name][i_last-1];
      for (let i = i_last; i < params[name].length; ++i)
        params[name][i] = val;
    }
  }
}

function updateParameterStructFromVector1(params, param_vec)
{
  let Nt = params.T_hist - tail_block_size;
  let num_blocks = Math.ceil(Nt / block_size);

  for (let k = 0; k < opt_array_names.length; ++k)
  {
    for (let i = 0; i < num_blocks; ++i)
    {
      for (let j = 0; j < block_size; ++j)
      {
        let ind = block_size*i + j;
        if (ind < Nt)
        {
          params[opt_array_names[k]][ind] = param_vec[k*num_blocks + i];
        }
      }
    }

    for (let i = Nt; i < params.T_hist; ++i)
      params[opt_array_names[k]][i] = params[opt_array_names[k]][Nt-1];
  }


  // let off = opt_array_names.length*num_blocks;
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

function getParameterBounds(params)
{
  return getParameterBounds1(params);
}

function getParameterBounds1(params)
{
  let Nt = params.T_hist - tail_block_size;
  let num_blocks = Math.ceil(Nt / block_size);

  let bounds = [];

  for (let i = 0; i < num_blocks; ++i)
    bounds.push({min: 0.0, max: 1.0, step: 1e-4}); //b1N

  for (let rep = 0; rep < opt_array_names.length-1; ++rep)
    for (let i = 0; i < num_blocks; ++i)
      bounds.push({min: 0.0, max: 1.0, step: 1e-4}); //ce, c0, c1

  // bounds.push({min: 1.0, max: 100.0, step: 1e-4}); //E0_0
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a0
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a10
  // bounds.push({min: 0.1, max: 1.0, step: 1e-4}); //a11
  // bounds.push({min: 0.0, max: 0.3, step: 1e-4}); //g0
  // bounds.push({min: 0.0, max: 0.3, step: 1e-4}); //g1
  // bounds.push({min: 0.0, max: 0.1, step: 1e-4}); //p1
  // bounds.push({min: 0.0, max: 0.3, step: 1e-4}); //g2
  // bounds.push({min: 0.0, max: 0.1, step: 1e-4}); //p2
  // bounds.push({min: 0.0, max: 0.3, step: 1e-4}); //g3
  // bounds.push({min: 0.0, max: 0.1, step: 1e-5}); //mu

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
  for (let i = 0; i < param_vec.length; ++i)
    param_vec[i] = bounds[i].min + Math.random()*(bounds[i].max - bounds[i].min);

  // for (let i = 0; i < param_vec.length; ++i)
  // {
  //   let step = 2*Math.random() - 1.0;
  //   param_vec[i] += 0.1*step*(bounds[i].max - bounds[i].min);
  //   param_vec[i] = Math.min(Math.max(param_vec[i], bounds[i].min), bounds[i].max);
  // }
}

function optimizeParameters()
{
  var t0 = performance.now();

  let T_pred_orig = sim_params.T_pred;
  let params = sim_params;
  params.T_pred = 0;

  let param_vec = getParameterVector(params);
  let param_bounds = getParameterBounds(params);

  updateParameterStructFromVector(params, param_vec);

  const n = param_vec.length;
  let dparam_vec = new Array(n).fill(0);
  let param_vec0 = new Array(n).fill(0);
  let param_vec_opt = new Array(n).fill(0);
  copyVector(param_vec, param_vec_opt);

  //Compute and store regularization matrix
  // params.reg_matrix = getDCTMatrix(params.T_hist);
  params.N_haar = Math.pow(2.0, Math.round(Math.log2(params.T_hist)));
  params.reg_matrix = getHaarMatrix(params.N_haar);

  //Calculate initial cost
  params.cost_scaling = 1.0;
  let cost = getOptCost(params, false);
  params.cost_scaling = 1.0/cost;
  // cost *= params.cost_scaling;
  cost = getOptCost(params);
  let cost_init = cost;

  let cost_rel_opt = 1.0;
  const cost_reduction_tol = 1e-4;
  const min_eta = 1e-6;

  for (let pass = 0; pass < 5; ++pass)
  {
    cost = getOptCost(params);

    let cost_prev = cost;
    let stalled_iter = 0;

    let iter = 0;
    for (iter = 0; iter < 100; ++iter)
    {
      let cost_grad = getOptCostGradient(params, param_bounds);

      trimUpdateVector(param_vec, cost_grad, param_bounds)

      //reset params struct, since Jacobian calc can change values slightly (by machine tolerances)
      updateParameterStructFromVector(params, param_vec);

      let eta0 = getLimitingStepSize(param_vec, cost_grad, param_bounds);
      console.log("  Iter " + iter + ": " + (cost/cost_init).toExponential(4) + ", eta0: " + eta0.toExponential(4));

      if (eta0 < 1e-13)
        break;

      for (let i = 0; i < n; ++i)
      {
        param_vec0[i] = param_vec[i];
        cost_grad[i] *= eta0;
      }

      let eta = 1.0;
      while (eta >= min_eta)
      {
        //Update param_vec based on (scaled) update vector
        for (let i = 0; i < n; ++i)
          param_vec[i] = param_vec0[i] - eta*cost_grad[i];

        //This checks for the validity of parameters, and may modify values
        updateParameterStructFromVector(params, param_vec);

        //Evaluate new cost
        let cost_new = getOptCost(params);
        // console.log("  " + eta.toExponential(4) + ", " + cost_new.toExponential(4));

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

  var t1 = performance.now()
  console.log("Optimization time: " + ((t1 - t0)/1000).toFixed(3) + "s.")

  updateParameters(true);
  return params;
}

function getOptCost(params, regularize = true)
{
  let sol_hist = predictModel(params);

  let wA = 1.0, wR = 1.0, wF = 1.0;


  let cost = 0.0, cost0 = 0.0;
  for (let i = 1; i < sol_hist.length; ++i)
  {
    //Error in number of active patients
    let num_active_pred = sol_hist[i].getNumActiveDiagnosed();
    // let num_active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];

    let num_recov_pred = sol_hist[i].getNumRecoveredDiagnosed();
    // let num_recov_true = data_real.categorized[i].y[1];

    let num_fatal_pred = sol_hist[i].getNumFatalDiagnosed();
    let num_fatal_true = data_real.categorized[i].y[2];

    let num_total_pred = num_active_pred + num_recov_pred + num_fatal_pred;
    let num_total_true = data_real.categorized[i].y[0];

    // num_active_pred = (num_active_pred > 0) ? Math.log(num_active_pred) : 0;
    // num_active_true = (num_active_true > 0) ? Math.log(num_active_true) : 0;
    // num_recov_pred = (num_recov_pred > 0) ? Math.log(num_recov_pred) : 0;
    // num_recov_true = (num_recov_true > 0) ? Math.log(num_recov_true) : 0;
    // num_fatal_pred = (num_fatal_pred > 0) ? Math.log(num_fatal_pred) : 0;
    // num_fatal_true = (num_fatal_true > 0) ? Math.log(num_fatal_true) : 0;
    // num_total_pred = (num_total_pred > 0) ? Math.log(num_total_pred) : 0;
    // num_total_true = (num_total_true > 0) ? Math.log(num_total_true) : 0;

    // let cA = (num_active_true == 0) ? 1 : (1.0/num_active_true);
    // let cR = (data_real.categorized[i].y[1] == 0) ? 1 : (1.0/data_real.categorized[i].y[1]);
    // let cF = (data_real.categorized[i].y[2] == 0) ? 1 : (1.0/data_real.categorized[i].y[2]);

    // let err_active = (num_active_pred - num_active_true);
    // let err_recov = (num_recov_pred - num_recov_true);
    let err_fatal = (num_fatal_pred - num_fatal_true);
    let err_total = (num_total_pred - num_total_true);

    // cost += wA*err_active*err_active + wR*err_recov*err_recov + wF*err_fatal*err_fatal;
    cost += wA*err_total*err_total + wF*err_fatal*err_fatal;
    cost0 += num_total_true*num_total_true;
  }

  cost *= params.cost_scaling;
  // cost /= cost0;

  //Add regularization terms
  if (regularize)
  {
    let N = params.N_haar; //sol_hist.length
    let reg_b1 = getL1Norm(getMatVecProd(params.reg_matrix, params.b1N, N));
    let reg_c0 = getL1Norm(getMatVecProd(params.reg_matrix, params.c0, N));
    let reg_c1 = getL1Norm(getMatVecProd(params.reg_matrix, params.c1, N));

    // let reg_b1 = getL1Norm(getDCT(params.b1N, sol_hist.length));
    // let reg_c0 = getL1Norm(getDCT(params.c0, sol_hist.length));
    // let reg_c1 = getL1Norm(getDCT(params.c1, sol_hist.length));
    const wgt = 1e-1;
    cost += wgt*reg_b1 + wgt*reg_c0 + wgt*reg_c1;
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

function trimUpdateVector(param_vec, dparam_vec, param_bounds)
{
  const n = dparam_vec.length;
  for (let i = 0; i < n; ++i)
  {
    if ( (param_vec[i] == param_bounds[i].min && -dparam_vec[i] < 0.0) ||
         (param_vec[i] == param_bounds[i].max && -dparam_vec[i] > 0.0) )
         dparam_vec[i] = 0.0;
  }
}

function getCurrentPredictionError()
{
  const active_ind = [0, 1, 4, 5, 6]; //E1d, I0d, I1d, I2d, I3d

  let res_sq = 0.0, res0_sq = 0.0;
  for (let i = 1; i < data_real.total.length; ++i)
  {
    let active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];
    let active_pred = 0;
    for (let j = 0; j < active_ind.length; ++j)
      active_pred += data_predicted.cat_diag[active_ind[j]][i].y;

    let err_a = active_pred - active_true;
    let err_r = data_predicted.cat_diag[2][i].y - data_real.categorized[i].y[1]; //recovered
    let err_d = data_predicted.cat_diag[3][i].y - data_real.categorized[i].y[2]; //fatal

    res_sq += err_a*err_a + err_r*err_r + err_d*err_d;
    res0_sq += active_true*active_true + data_real.categorized[i].y[1]*data_real.categorized[i].y[1] + data_real.categorized[i].y[2]*data_real.categorized[i].y[2];
  }
  return Math.sqrt(res_sq/res0_sq);
}

function getL1Norm(vec)
{
  let norm = 0.0;
  for (let i = 0; i < vec.length; ++i)
    norm += Math.abs(vec[i]);
  return norm;
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

function getDCTMatrix(N)
{
  const scale_row0 = 1.0 / Math.sqrt(N);
  const scale_rowk = Math.sqrt(2.0) / Math.sqrt(N);

  let A = new Array(N*N).fill(0);

  //A(0,:)
  for (let n = 0; n < N; ++n)
    A[n] = scale_row0; //row k = 0

  //A(1:N,:)
  for (let k = 1; k < N; ++k)
    for (let n = 0; n < N; ++n)
      A[k*N + n] = Math.cos(Math.PI*(n + 0.5)*k / N) * scale_rowk;

  return A;
}

//Computes the discrete cosine transform of the vector x
function getDCT(x, N = x.length)
{
  const sqrtN = Math.sqrt(N);
  const sqrt2 = Math.sqrt(2.0);

  let X = new Array(N).fill(0);

  for (let k = 0; k < N; ++k)
  {
    let sum = 0.0;
    for (let n = 0; n < N; ++n)
      sum += x[n]*Math.cos(Math.PI*(n + 0.5)*k / N);

    if (k == 0)
      X[k] = sum / sqrtN;
    else
      X[k] = sum / sqrtN * sqrt2;
  }
  return X;
}

function getHaarMatrix(N)
{
  let num_levels = Math.floor(Math.log2(N));
  if (Math.pow(2,num_levels) != N)
    throw("N is not a power of 2!");

  let A = new Array(N*N).fill(0);
  const inv_sqrtN = 1.0 / Math.sqrt(N);

  for (let j = 0; j < N; ++j)
    A[j] = inv_sqrtN;

  for (let k = 1; k < N; ++k)
  {
    let p = Math.floor(Math.log(k) / Math.log(2.0));
    let k1 = Math.pow(2.0, p);
    let k2 = k1 * 2.0;
    let q = k - k1;
    let t1 = N/k1;
    let t2 = N/k2;
    let tmp = Math.pow(2.0, p/2.0) * inv_sqrtN;
    for (let i = 0; i < t2; ++i)
    {
      A[N*k + q*t1 + i] = tmp;
      A[N*k + q*t1 + i + t2] = -tmp;
    }
  }

  return A;
}

function getMatVecProd(A, x, n = x.length)
{
  const m = A.length / n;
  let y = new Array(m).fill(0);

  for (let i = 0; i < m; ++i)
    for (let j = 0; j < n; ++j)
      y[i] += A[i*n + j] * x[j];
  return y;
}
