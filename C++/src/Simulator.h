#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>

class ModelParams;
class Population;

std::vector<Population> predictModel(const ModelParams& params, const Population& pop_init);

#endif
