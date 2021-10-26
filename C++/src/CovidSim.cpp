#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <cmath>

#include "models/Simulator.h"
#include "models/Population.h"
#include "models/ModelParams.h"

int
main(int argc, char *argv[])
{
  if (argc == 1 || argc > 3)
    throwError("Invalid arguments! Argument format: country_name [population = 1e7]");

  std::string country;
  int N_population = 1e7;

  if (argc >= 2)
    country = argv[1];
  if (argc >= 3)
    N_population = std::stod(argv[2]);

  std::cout << "Country name: " << country << std::endl;
  std::cout << "Population: " << N_population << std::endl;

  std::string filepath = "csv_data/" + country + "_params.txt";
  std::ifstream in(filepath);
  if (!in)
    throwError("Cannot open file - " + filepath);

  std::vector<double> param_vec;
  double val;
  while (in >> val)
    param_vec.push_back(val);
  in.close();

	Population<double> pop_init(N_population, 5);

	std::string folder_path = "results";
	mkdir(folder_path.c_str(), 0777);

  std::string filepath_predictions = folder_path + "/" + country + ".txt";
  std::ofstream file_predictions(filepath_predictions);
  if (!file_predictions.good())
    throwError("Cannot open file to write - " + filepath_predictions);

  const int nt_hist = (param_vec.size() - 11) / 5;
  if (nt_hist <= 0)
    throwError("nt_hist <= 0!");

  ModelParams<double> params(nt_hist, 0);
  copyFullVector2Param<double>(param_vec, params);

  auto t0 = std::chrono::high_resolution_clock::now();

  auto prediction = predictModel(params, pop_init);

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << "ms" << std::endl;

  file_predictions << N_population << std::endl;
  for (std::size_t i = 0; i < prediction.size(); ++i)
  {
    file_predictions << (int) std::round(prediction[i].getNumReported()) << " "
                     << (int) std::round(prediction[i].getNumRecoveredReported()) << " "
                     << (int) std::round(prediction[i].getNumFatalReported()) << std::endl;
  }
  std::cout << "Wrote prediction to " << filepath_predictions << std::endl;
  file_predictions.close();

	return 0;
}

