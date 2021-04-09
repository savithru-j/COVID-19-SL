#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <iomanip>

#include "Simulator.h"
#include "OptimizerFull.h"

int
main(int argc, char *argv[])
{
  if (argc == 1 || argc > 9)
    throwError("Invalid arguments! Argument format: country_name "
               "[w_conf = 1.0] [w_recov = 1.0] [w_fatal = 1.0] [w_reg = 0.01] "
               "[max_iter_per_pass = 1000] [max_passes = 1] [seed = 1]");

  std::string country;
  double weight_conf = 1.0, weight_recov = 1.0, weight_fatal = 1.0, weight_reg = 0.01;
  int max_iter_per_pass = 1000, max_passes = 1;
  int seed = 1;

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
  if (argc >= 9)
    seed = std::stoi(argv[8]);

  std::cout << "Country name: " << country << std::endl;
  std::cout << "Optimization weights: " << std::endl;
  std::cout << " w_confirmed: " << weight_conf << std::endl;
  std::cout << " w_recovered: " << weight_recov << std::endl;
  std::cout << " w_fatal: " << weight_fatal << std::endl;
  std::cout << " w_regularization: " << weight_reg << std::endl;
  std::cout << "Max iterations per pass: " << max_iter_per_pass << std::endl;
  std::cout << "Max passes: " << max_passes << std::endl;
  std::cout << "Random seed: " << seed << std::endl;
  std::cout << std::endl;

	ObservedPopulation pop_observed("csv_data/" + country + ".txt");

	Population pop_init(pop_observed.N, 5, 1);

	std::string folder_path = "results";
	mkdir(folder_path.c_str(), 0777);

  std::string suffix = "_seed" + std::to_string(seed);
	std::string filepath_opt_params = folder_path + "/" + country + "_params" + suffix + ".txt";
  std::ofstream file_opt_params(filepath_opt_params);
  if (!file_opt_params.good())
    throwError("Cannot open file to write - " + filepath_opt_params);

  std::string filepath_predictions = folder_path + "/" + country + "_prediction" + suffix + ".txt";
  std::ofstream file_predictions(filepath_predictions);
  if (!file_predictions.good())
    throwError("Cannot open file to write - " + filepath_predictions);

  auto t0 = std::chrono::high_resolution_clock::now();

  std::string filepath_quarantine_input = "csv_data/" + country + "_external.txt";
  std::cout << "Checking for external quarantine input data at " << filepath_quarantine_input << std::endl;
  auto pop_quarantine_input = getQuarantineInputVector(filepath_quarantine_input);
  std::cout << "External quarantine input vector length: " << pop_quarantine_input.size() << std::endl << std::endl;

  OptimizerFull opt(pop_observed, pop_init, pop_quarantine_input,
                    weight_conf, weight_recov, weight_fatal, weight_reg,
                    max_iter_per_pass, max_passes, seed);

#if 0
  std::string filepath = "csv_data/" + country + "_params.txt";
  std::ifstream in(filepath);
  std::vector<double> param_vec;
  double val;
  while (in >> val)
    param_vec.push_back(val);
  in.close();
  opt.param_vec = param_vec;
  Optimizer::copyVector2Param(opt.param_vec, opt.params);
#endif

  opt.optimizeParametersNLOPT();

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;

  const int nt_hist = pop_observed.getNumDays();
  std::array<std::vector<Population>, OptimizerFull::NUM_RESULTS> predictions;
  std::array<Vector, OptimizerFull::NUM_RESULTS> param_vecs_full;

  for (std::size_t i = 0; i < predictions.size(); ++i)
  {
    ModelParams params(nt_hist, opt.t_buffer);
    opt.copyVector2Param(opt.optimal_param_vec[i], params);
    predictions[i] = predictModel(params, pop_init);

    param_vecs_full[i].resize(5*(params.nt_hist + opt.t_buffer) + 10);
    copyParam2FullVector(params, param_vecs_full[i]);
  }

  file_opt_params << std::scientific << std::setprecision(6);
  for (int i = 0; i < param_vecs_full[0].m(); ++i)
  {
    for (std::size_t j = 0; j < param_vecs_full.size()-1; ++j)
      file_opt_params << param_vecs_full[j][i] << ", ";
    file_opt_params << param_vecs_full.back()[i] << std::endl;
  }
  file_opt_params << opt.cost_min << std::endl;
  std::cout << "Wrote optimal parameters to " << filepath_opt_params << std::endl;


  file_predictions << std::scientific << std::setprecision(6);
  for (std::size_t i = 0; i < predictions[0].size(); ++i)
  {
    for (std::size_t j = 0; j < predictions.size(); ++j)
    file_predictions << predictions[j][i].getNumReported() << ", "
                     << predictions[j][i].getNumRecoveredReported() << ", "
                     << predictions[j][i].getNumFatalReported() << ", "
                     << predictions[j][i].getNumInfectedUnreported() << ", "
                     << predictions[j][i].getNumRecoveredUnreported() << ", "
                     << predictions[j][i].getNumFatalUnreported() << ", ";
    file_predictions << std::endl;
  }
  std::cout << "Wrote predictions to " << filepath_predictions << std::endl;

  file_opt_params.close();
  file_predictions.close();

	return 0;
}

