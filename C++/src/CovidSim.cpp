#include <iostream>
#include <chrono>

#include "Population.h"
#include "LinearAlgebra.h"

std::vector<Population> predictModel(const ModelParams& params, const Population& pop_init);

int main()
{

  ModelParams params(75, 14);

  Population pop0(21323733, 5, 1);

  auto t0 = std::chrono::high_resolution_clock::now();

  std::vector<Population> pop_hist;
  for (int i = 0; i < 1e1; i++)
    pop_hist = predictModel(params, pop0);

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count()/1000.0 << "s" << std::endl;

  std::cout << pop_hist.size() << std::endl;
	std::cout << pop_hist.back() << std::endl;

	Matrix A(5,6);
	Vector x = std::vector<double>({1,2,3,4,5,6});
	auto y = A*x;

	ObservedPopulation("csv_data/srilanka.txt");

	return 0;
}


std::vector<Population> predictModel(const ModelParams& params, const Population& pop_init)
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
