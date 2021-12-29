#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

template<class T>
class ModelParams;

template<class T>
class ModelParams2Layer;

template<class T>
class ModelParams4Layer;

template<class T>
class Population;

template<class T>
class Population2Layer;

template<class T>
class Population4Layer;

template<template<class> class ParamsT, template<class> class PopulationT, class T>
std::vector<PopulationT<T>>
predictModel(const ParamsT<T>& params, const PopulationT<double>& pop_init)
{
  std::vector<PopulationT<T>> population_hist;

  const int nt = params.nt_hist + params.nt_pred;
  const int nt_sub = 1.0/params.dt;

  PopulationT<T> pop = pop_init;

  population_hist.reserve(nt);
  population_hist.push_back(pop);

  for (int t = 0; t < nt-1; t++)
  {
    for (int j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop.evolve(params, t);
    pop.report(params, t); //report once daily
    pop.vaccinate(params, t); //remove vaccinated individuals daily

    population_hist.push_back(pop); //save solution daily
  }
  return population_hist;
}

template<class T>
T calcEffectiveReproductionRatio(const ModelParams<T>& params_orig,
                                 const std::vector<Population<T>>& pop_hist_orig, const int t);

template<class T>
T calcEffectiveReproductionRatio(const ModelParams2Layer<T>& params_orig,
                                 const std::vector<Population2Layer<T>>& pop_hist_orig, const int t);

template<class T>
T calcEffectiveReproductionRatio(const ModelParams4Layer<T>& params_orig,
                                 const std::vector<Population4Layer<T>>& pop_hist_orig, const int t);

#endif
