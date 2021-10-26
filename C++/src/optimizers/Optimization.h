#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <nlopt.hpp>

struct ParamBound {

  ParamBound() = default;
  ParamBound(double min_, double max_, double step_) : min(min_), max(max_), step(step_) {}

  double min = 0.0, max = 1.0, step = 1e-4;
};

std::string getNLOPTResultDescription(nlopt::result resultcode);

#endif
