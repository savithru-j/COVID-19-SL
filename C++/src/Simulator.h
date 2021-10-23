#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

template<class T>
class ModelParams;

template<class T>
class Population;

template<class T>
std::vector<Population<T>>
predictModel(const ModelParams<T>& params, const Population<T>& pop_init);

template<class T>
T calcEffectiveReproductionRatio(const ModelParams<T>& params_orig,
                                 const std::vector<Population<T>>& pop_hist_orig, const int t);

#endif
