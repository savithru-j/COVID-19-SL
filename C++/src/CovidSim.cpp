#include <iostream>
#include <chrono>

#include "CovidSim.h"

int
main()
{
  ModelParams params(75, 14);

  Population pop0(21323733, 5, 1);

  auto t0 = std::chrono::high_resolution_clock::now();

  std::vector<Population> pop_hist;
  pop_hist = predictModel(params, pop0);

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;

  for (int i = 0; i < 10; i++)
    std::cout << uniformRand() << std::endl;

	ObservedPopulation pop_SL("csv_data/srilanka.txt");

	Vector x = std::vector<double>{0.4, 5, 2.7, -0.3, 3.6};
	Matrix A = getHaarMatrix(x.size());
	Vector y = A*x;

	getOptimalParameters(pop_SL);

	return 0;
}


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

    population_hist.push_back(pop); //save solution daily
  }

  return population_hist;
}

ModelParams
getOptimalParameters(const ObservedPopulation& observed_pop)
{
  const int T_end_buffer = 5; //no. of non-optimized days at end

  const int nt = observed_pop.getNumDays();
  const int nt_opt = nt - T_end_buffer; //number of days to optimize parameters for

  const auto param_bounds = getParameterBounds(nt_opt);
  Vector param_vec(param_bounds.size());

  double jump_energy = 1.0;
  randomizeParameterVector(param_bounds, jump_energy, param_vec);

  ModelParams params(nt_opt, T_end_buffer);
  copyVector2Param(param_vec, params);

  Population pop_init(observed_pop.N, 5, 1); //TODO: set initial E0 and Rd



  double cost = getOptCost(params,param_vec, pop_init, observed_pop);


  return params;
}

void copyParam2Vector(const ModelParams& params, Vector& v)
{
  const int nt = params.nt_hist;
  const int m = 3*nt + 10;
  if (v.m() != m)
    v.resize(m);

  for (int i = 0; i < nt; ++i)
  {
    v[       i] = params.betaN[i];
    v[  nt + i] = params.ce[i];
    v[2*nt + i] = params.c1[i];
  }
  const int off = 3*nt;
  v[off  ] = params.T_incub0;
  v[off+1] = params.T_incub1;
  v[off+2] = params.T_asympt;
  v[off+3] = params.T_mild;
  v[off+4] = params.T_severe;
  v[off+5] = params.T_icu;
  v[off+6] = params.f;
  v[off+7] = params.frac_recover_I1;
  v[off+8] = params.frac_recover_I2;
  v[off+9] = params.CFR;
}

void copyVector2Param(const Vector& v, ModelParams& params)
{
  const int nt = params.nt_hist;
  const int m = 3*nt + 10;
  if (v.m() != m)
    throwError("copyVector2Param - inconsistent dimensions!");

  for (int i = 0; i < nt; ++i)
  {
    params.betaN[i] = v[       i];
    params.ce[i]    = v[  nt + i];
    params.c1[i]    = v[2*nt + i];
  }

  const int N = params.nt_hist + params.nt_pred;
  for (int i = nt; i < N; ++i) //Copy last "history" value to prediction section
  {
    params.betaN[i] = v[  nt-1];
    params.ce[i]    = v[2*nt-1];
    params.c1[i]    = v[3*nt-1];
  }

  const int off = 3*nt;
  params.T_incub0        = v[off  ];
  params.T_incub1        = v[off+1];
  params.T_asympt        = v[off+2];
  params.T_mild          = v[off+3];
  params.T_severe        = v[off+4];
  params.T_icu           = v[off+5];
  params.f               = v[off+6];
  params.frac_recover_I1 = v[off+7];
  params.frac_recover_I2 = v[off+8];
  params.CFR             = v[off+9];
}

std::vector<ParamBound>
getParameterBounds(int nt)
{
  const int m = 3*nt + 10;
  std::vector<ParamBound> bounds(m);

  for (int i = 0; i < nt; ++i)
  {
    bounds[       i] = ParamBound(0, 1, 1e-4); //betaN
    bounds[  nt + i] = ParamBound(0, 1, 1e-4); //ce
    bounds[2*nt + i] = ParamBound(0, 1, 1e-4); //c1
  }
  const int off = 3*nt;
  bounds[off  ] = ParamBound(1, 10, 1e-4); //T_incub0
  bounds[off+1] = ParamBound(1, 10, 1e-4); //T_incub1
  bounds[off+2] = ParamBound(1, 10, 1e-4); //T_asympt
  bounds[off+3] = ParamBound(1, 10, 1e-4); //T_mild
  bounds[off+4] = ParamBound(1, 10, 1e-4); //T_severe
  bounds[off+5] = ParamBound(1, 10, 1e-4); //T_icu
  bounds[off+6] = ParamBound(0, 1, 1e-4); //f
  bounds[off+7] = ParamBound(0, 1, 1e-4); //frac_recover_I1
  bounds[off+8] = ParamBound(0, 1, 1e-4); //frac_recover_I2
  bounds[off+9] = ParamBound(0, 0.02, 1e-5); //CFR

  return bounds;
}

double uniformRand(double min, double max)
{
//  static std::random_device rd;  //Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(2.0); //Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> uniform_rand(min, max);
  return uniform_rand(gen);
}

void
randomizeParameterVector(const std::vector<ParamBound>& bounds, const double energy, Vector& param_vec)
{
  for (int i = 0; i < param_vec.m(); ++i)
  {
    const double delta = energy * (bounds[i].max - bounds[i].min) * uniformRand(-1, 1);
    param_vec[i] = std::min(std::max(param_vec[i] + delta, bounds[i].min), bounds[i].max);
  }
}

double getOptCost(const ModelParams& params, const Vector& param_vec, const Population& pop_init,
                  const ObservedPopulation& observed_pop, const Matrix& reg_matrix,
                  const double scaling)
{
  auto pop_hist = predictModel(params, pop_init);

  const double wC = 1.0, wR = 1.0, wF = 1.0;
  const int nt = params.nt_hist + params.nt_pred;

  double cost = 0.0;
  for (int i = 1; i < nt; ++i)
  {
    double err_total = (pop_hist[i].getNumReported() - observed_pop.confirmed[i]);
    double err_recov = (pop_hist[i].getNumRecoveredReported() - observed_pop.recovered[i]);
    double err_fatal = (pop_hist[i].getNumFatalReported() - observed_pop.deaths[i]);

    cost += wC*err_total*err_total + wR*err_recov*err_recov + wF*err_fatal*err_fatal;
  }

  cost *= scaling;

  //Add regularization terms
  if (reg_matrix.m() > 0 && reg_matrix.n() > 0)
  {
    double reg_betaN = 0.0;
    double reg_ce = 0.0;
    double reg_c1 = 0.0;

    for (int i = 0; i < reg_matrix.m(); ++i)
    {
      double coeff_betaN = 0.0;
      double coeff_ce = 0.0;
      double coeff_c1 = 0.0;

      for (int j = 0; j < reg_matrix.n(); ++j)
      {
        if (j < params.betaN.m())
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN(j);
          coeff_ce += reg_matrix(i,j) * params.ce(j);
          coeff_c1 += reg_matrix(i,j) * params.c1(j);
        }
        else
        {
          coeff_betaN += reg_matrix(i,j) * params.betaN.back();
          coeff_ce += reg_matrix(i,j) * params.ce.back();
          coeff_c1 += reg_matrix(i,j) * params.c1.back();
        }
      }

      //Take L1-norm of coefficient vector
      reg_betaN += std::abs(coeff_betaN);
      reg_ce += std::abs(coeff_ce);
      reg_c1 += std::abs(coeff_c1);
    }

    const double wgt = 1e-1;
    cost += wgt*reg_betaN + wgt*reg_ce + wgt*reg_c1;
  }

//  if (isNaN(cost))
//  {
//    console.log("getOptimizationCost: found NaN");
//    return {};
//  }
  return cost;
}

Matrix getHaarMatrix(int m)
{
  const int num_levels = std::ceil(std::log2(m));
  const int n = std::pow(2, num_levels);
//  if (std::pow(2,num_levels) != n)
//    throw("N is not a power of 2!");

  Matrix A(n, n, 0.0);
  const double inv_sqrtn = 1.0 / std::sqrt(n);

  for (int j = 0; j < n; ++j)
    A(0, j) = inv_sqrtn;

  for (int k = 1; k < n; ++k)
  {
    int p = std::floor(std::log(k) / std::log(2.0));
    double k1 = std::pow(2.0, p);
    double k2 = k1 * 2.0;
    double q = k - k1;
    double t1 = n/k1;
    double t2 = n/k2;
    double tmp = std::pow(2.0, p/2.0) * inv_sqrtn;
    for (int i = 0; i < t2; ++i)
    {
      A(k, q*t1 + i) = tmp;
      A(k, + q*t1 + i + t2) = -tmp;
    }
  }

  return A;
}
