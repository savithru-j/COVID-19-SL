#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <iomanip>

#include <nlopt.hpp>

#include "Simulator.h"
#include "Optimization.h"

int
main()
{
  std::string country = "srilanka64";

	ObservedPopulation pop_observed("csv_data/" + country + ".txt");

	Population pop_init(pop_observed.N, 5, 1);

	std::string folder_path = "results";
	mkdir(folder_path.c_str(), 0777);

	std::string filepath_opt_params = folder_path + "/" + country + "_params.txt";
  std::ofstream file_opt_params(filepath_opt_params);
  if (!file_opt_params.good())
    throwError("Cannot open file to write - " + filepath_opt_params);

  std::string filepath_predictions = folder_path + "/" + country + "_prediction.txt";
  std::ofstream file_predictions(filepath_predictions);
  if (!file_predictions.good())
    throwError("Cannot open file to write - " + filepath_predictions);

  auto t0 = std::chrono::high_resolution_clock::now();

//  auto optinfo2 = optimizeParametersNLOPT(pop_observed, pop_init);
//  exit(EXIT_FAILURE);

  double reg_weight = 0.01;
  Optimizer opt(pop_observed, pop_init, reg_weight);
//  opt.optimizeParameters();
  opt.optimizeParametersNLOPT();

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;


  file_opt_params << std::scientific << std::setprecision(6);
  for (int i = 0; i < opt.optimal_param_vec[0].m(); ++i)
  {
    for (std::size_t j = 0; j < opt.optimal_param_vec.size()-1; ++j)
      file_opt_params << opt.optimal_param_vec[j][i] << ", ";
    file_opt_params << opt.optimal_param_vec.back()[i] << std::endl;
  }

  const int nt_hist = pop_observed.getNumDays();

  std::array<std::vector<Population>, Optimizer::NUM_RESULTS> predictions;
  for (std::size_t i = 0; i < predictions.size(); ++i)
  {
    ModelParams params(nt_hist - opt.t_buffer, opt.t_buffer);
    copyVector2Param(opt.optimal_param_vec[i], params);
    predictions[i] = predictModel(params, pop_init);
  }

  file_predictions << std::scientific << std::setprecision(6);
  for (std::size_t i = 0; i < predictions[0].size(); ++i)
  {
    for (std::size_t j = 0; j < predictions.size(); ++j)
    file_predictions << predictions[j][i].getNumReported() << ", "
                     << predictions[j][i].getNumRecoveredReported() << ", "
                     << predictions[j][i].getNumFatalReported() << ", ";
    file_predictions << std::endl;
  }

  file_opt_params.close();
  file_predictions.close();

//  ModelParams params(nt_hist - 5, 0);
//  copyVector2Param(optinfo.param_vec_opt[0], params);
//  Vector x(128);
//  for (int i = 0; i < 128; ++i)
//  {
//    if (i < params.betaN.m())
//      x[i] = params.betaN[i];
//    else
//      x[i] = params.betaN.back();
//  }
//
//  Vector coeff = optinfo.reg_matrix*x;
//
//  std::cout << "coeff: " << coeff << std::endl;

	return 0;
}

