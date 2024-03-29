#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <iomanip>

#include "models/Simulator.h"
#include "optimizers/OptimizerPiecewise2Layer.h"

int
main(int argc, char *argv[])
{
  if (argc == 1 || argc > 10)
    throwError("Invalid arguments! Argument format: country_name "
               "[w_conf = 1.0] [w_recov = 1.0] [w_fatal = 1.0] [interval_size = 5] [linear_basis = 0] "
               "[max_iter_per_pass = 1000] [max_passes = 1] [seed = 1]");

  std::string country;
  double weight_conf = 1.0, weight_recov = 1.0, weight_fatal = 1.0;
  int interval_size = 5, linear_basis = 0, max_iter_per_pass = 1000, max_passes = 1;
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
    interval_size = std::stoi(argv[5]);
  if (argc >= 7)
    linear_basis = std::stoi(argv[6]);
  if (argc >= 8)
    max_iter_per_pass = std::stoi(argv[7]);
  if (argc >= 9)
    max_passes = std::stoi(argv[8]);
  if (argc >= 10)
    seed = std::stoi(argv[9]);

  std::cout << "Country name: " << country << std::endl;
  std::cout << "Optimization weights: " << std::endl;
  std::cout << " w_confirmed: " << weight_conf << std::endl;
  std::cout << " w_recovered: " << weight_recov << std::endl;
  std::cout << " w_fatal: " << weight_fatal << std::endl;
  std::cout << "Interval size: " << interval_size << " days" << std::endl;
  std::cout << "Linear basis? " << linear_basis << std::endl;
  std::cout << "Max iterations per pass: " << max_iter_per_pass << std::endl;
  std::cout << "Max passes: " << max_passes << std::endl;
  std::cout << "Random seed: " << seed << std::endl;
  std::cout << std::endl;

  std::string input_data_folder = "csv_data/";
	ObservedPopulation pop_observed(input_data_folder + country + ".txt");

	std::string folder_path = "results";
	mkdir(folder_path.c_str(), 0777);

	int T_smooth = 7;
  auto confirmed_vec = pop_observed.confirmed; //save original data before smoothing
  auto death_vec = pop_observed.deaths;
  pop_observed.smooth(T_smooth);

  Population2Layer<double> pop_init(pop_observed.N, 100,
                                    pop_observed.confirmed[0] - pop_observed.recovered[0] - pop_observed.deaths[0],
                                    pop_observed.recovered[0],
                                    pop_observed.deaths[0]);

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

  std::string filepath_vaccine_data = input_data_folder + country + "_vaccinations.txt";
  std::cout << "Checking for vaccination data at " << filepath_vaccine_data << std::endl;
  auto pop_vaccine_data = readDataVectorFromFile(filepath_vaccine_data, T_smooth);
  std::cout << "Vaccination data vector length: " << pop_vaccine_data.size() << std::endl << std::endl;

  if (pop_vaccine_data.size() > 0 && (int) pop_vaccine_data.size() != pop_observed.getNumDays())
    throwError("Length of vaccine data does not match observed population data!");

  std::string filepath_testing_data = input_data_folder + country + "_testing.txt";
  std::cout << "Checking for testing data at " << filepath_testing_data << std::endl;
  auto pop_testing_data = readDataVectorFromFile(filepath_testing_data, T_smooth);
  std::cout << "Testing data vector length: " << pop_testing_data.size() << std::endl << std::endl;

  if (pop_testing_data.size() > 0 && (int) pop_testing_data.size() != pop_observed.getNumDays())
    throwError("Length of testing data does not match observed population data!");

  std::ofstream file_popsmooth(folder_path + "/" + country + "_popsmooth.txt");
  for (std::size_t i = 0; i < pop_observed.confirmed.size(); ++i)
  {
    file_popsmooth << confirmed_vec[i] << ", " << death_vec[i] << ", "
                   << pop_observed.confirmed[i] << ", " << pop_observed.deaths[i];
    if (pop_vaccine_data.size() > 0)
      file_popsmooth << ", " << (int) pop_vaccine_data[i];
    if (pop_testing_data.size() > 0)
      file_popsmooth << ", " << (int) pop_testing_data[i];
    file_popsmooth << std::endl;
  }
  file_popsmooth.close();

  OptimizerPiecewise2Layer opt(pop_observed, pop_init, pop_testing_data, pop_vaccine_data,
                               interval_size, linear_basis, weight_conf, weight_recov, weight_fatal,
                               max_iter_per_pass, max_passes, seed);

#if 0 //Initialize to true params
  std::string filepath = "csv_data/" + country + "_params.txt";
  std::ifstream in(filepath);
  std::vector<double> param_vec;
  double val;
  while (in >> val)
    param_vec.push_back(val);
  in.close();
  Optimizer::copyVector2Param(param_vec, opt.params);
  opt.copyParam2Vector(opt.params, opt.param_vec);
#endif

  opt.optimizeParametersNLOPT();

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;

  const int nt_hist = pop_observed.getNumDays();
  std::vector<std::vector<Population2Layer<double>>> predictions(opt.cost_min.size());
  std::vector<Vector<double>> param_vecs_full(opt.cost_min.size());

  std::vector<Vector<double>> Reff_vector(predictions.size(), Vector<double>(nt_hist));

  for (std::size_t i = 0; i < predictions.size(); ++i)
  {
    ModelParams2Layer<double> params(nt_hist, 0, pop_testing_data, pop_vaccine_data);
    opt.copyVector2Param(opt.optimal_param_vec[i], params);
    predictions[i] = predictModel(params, pop_init);

    param_vecs_full[i].resize(2*params.nt_hist + 6);
    copyParam2FullVector(params, param_vecs_full[i]);

    //Calculate R-eff
    for (std::size_t t = 0; t < predictions[0].size(); ++t)
      Reff_vector[i][t] = calcEffectiveReproductionRatio(params, predictions[i], t);
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

  //Write header
  for (std::size_t j = 0; j < predictions.size(); ++j)
  {
    file_predictions << "TotalR" << j
                     << ", TotalU" << j
                     << ", FatalR" << j
                     << ", InfectedR" << j
                     << ", InfectedU" << j
                     << ", TotalVac" << j
                     << ", Reff" << j;
    if (j < predictions.size()-1)
      file_predictions << ", ";
  }
  file_predictions << std::endl;

  file_predictions << std::scientific << std::setprecision(10);
  for (std::size_t t = 0; t < predictions[0].size(); ++t)
  {
    for (std::size_t isol = 0; isol < predictions.size(); ++isol)
    {
      file_predictions << predictions[isol][t].getNumReported() << ", "
                       << predictions[isol][t].getNumUnreported() << ", "
                       << predictions[isol][t].getNumFatalReported() << ", "
                       << predictions[isol][t].getNumInfectedReported() << ", "
                       << predictions[isol][t].getNumInfectedUnreported() << ", "
                       << predictions[isol][t].getNumVaccinated() << ", "
                       << Reff_vector[isol][t];
      if (isol < predictions.size()-1)
        file_predictions << ", ";
    }
    file_predictions << std::endl;
  }
  std::cout << "Wrote predictions to " << filepath_predictions << std::endl;

  file_opt_params.close();
  file_predictions.close();

	return 0;
}

