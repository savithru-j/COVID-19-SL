#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <iomanip>

#include <nlopt.hpp>

#include "Simulator.h"
#include "Optimization.h"

int
main(int argc, char *argv[])
{
  if (argc == 1 || argc > 8)
    throwError("Invalid arguments! Argument format: country_name "
               "[w_conf = 1.0] [w_recov = 1.0] [w_fatal = 1.0] [w_reg = 0.01] [max_iter_per_pass = 1000] [max_passes = 1]");

  std::string country;
  double weight_conf = 1.0, weight_recov = 1.0, weight_fatal = 1.0, weight_reg = 0.01;
  int max_iter_per_pass = 1000, max_passes = 1;

  if (argc >= 2)
    country = argv[1];
  if (argc >= 3)
    weight_conf = std::stod(argv[2]);
  if (argc >= 4)
    weight_recov = std::stod(argv[3]);
  if (argc >= 5)
    weight_fatal = std::stod(argv[4]);
  if (argc >= 6)
    weight_reg = std::stod(argv[5]);
  if (argc >= 7)
    max_iter_per_pass = std::stoi(argv[6]);
  if (argc >= 8)
    max_passes = std::stoi(argv[7]);

  std::cout << "Country name: " << country << std::endl;
  std::cout << "Optimization weights: " << std::endl;
  std::cout << " w_confirmed: " << weight_conf << std::endl;
  std::cout << " w_recovered: " << weight_recov << std::endl;
  std::cout << " w_fatal: " << weight_fatal << std::endl;
  std::cout << " w_regularization: " << weight_reg << std::endl;
  std::cout << "Max iterations per pass: " << max_iter_per_pass << std::endl;
  std::cout << "Max passes: " << max_passes << std::endl;
  std::cout << std::endl;

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

  Optimizer opt(pop_observed, pop_init, weight_conf, weight_recov, weight_fatal, weight_reg,
                max_iter_per_pass, max_passes);
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
  std::cout << "Wrote optimal parameters to " << filepath_opt_params << std::endl;

  const int nt_hist = pop_observed.getNumDays();

  std::array<std::vector<Population>, Optimizer::NUM_RESULTS> predictions;
  for (std::size_t i = 0; i < predictions.size(); ++i)
  {
    ModelParams params(nt_hist - opt.t_buffer, opt.t_buffer + 7);
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
  std::cout << "Wrote predictions to " << filepath_predictions << std::endl;

  file_opt_params.close();
  file_predictions.close();

	return 0;
}

