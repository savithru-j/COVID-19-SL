#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

template<class T>
class ModelParams;

template<class T>
class ModelParams4Layer;

template<class T>
class Population;

template<class T>
class Population4Layer;

template<class T>
std::vector<Population<T>>
predictModel(const ModelParams<T>& params, const Population<double>& pop_init);

template<class T>
std::vector<Population4Layer<T>>
predictModel(const ModelParams4Layer<T>& params, const Population4Layer<double>& pop_init);

template<class T>
T calcEffectiveReproductionRatio(const ModelParams<T>& params_orig,
                                 const std::vector<Population<T>>& pop_hist_orig, const int t);

template<class T>
T calcEffectiveReproductionRatio(const ModelParams4Layer<T>& params_orig,
                                 const std::vector<Population4Layer<T>>& pop_hist_orig, const int t);

#endif
