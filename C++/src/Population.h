#ifndef POPULATION_H
#define POPULATION_H

#include <array>
#include <fstream>

#include "ModelParams.h"
#include "ErrorHandler.h"

class Population
{
public:

  using Array2 = std::array<double,2>;

  Population(double N_, double E0_0 = 0.0, double Rd_0 = 0.0)
    : N(N_), S(N - E0_0 - Rd_0), E0(E0_0), E1({0,0}),
      I0({0,0}), I1({0,0}), I2({0,0}), I3({0,0}), R({0, Rd_0}), D({0,0}) {}

  void evolve(const ModelParams& params, int t);
  void report(const ModelParams& params, int t);

  inline double getNumReported() const { return E1[1] + I0[1] + I1[1] + I2[1] + I3[1] + R[1] + D[1]; }
  inline double getNumActiveReported() const { return E1[1] + I0[1] + I1[1] + I2[1] + I3[1] + R[1]; }
  inline double getNumRecoveredReported() const { return R[1]; }
  inline double getNumFatalReported() const { return D[1]; }

  inline double getNumInfectedUnreported() const { return I0[0] + I1[0] + I2[0] + I3[0]; }
  inline double getNumRecoveredUnreported() const { return R[0]; }
  inline double getNumFatalUnreported() const { return D[0]; }

  double N = 0;
  double S = 0;
  double E0 = 0;
  Array2 E1;
  Array2 I0;
  Array2 I1;
  Array2 I2;
  Array2 I3;
  Array2 R;
  Array2 D;
  double dS_exit_Reff = 0;
};

std::ostream& operator<<(std::ostream& os, const Population& pop);


class ObservedPopulation
{
public:

  ObservedPopulation(const std::string& filepath)
  {
      std::ifstream in(filepath);
      if (!in)
        throwError("Cannot open file - " + filepath);

      in >> N; //read in population
      int c;
      while (in >> c)
      {
        int r, d;
        in >> r >> d;
        confirmed.push_back(c);
        recovered.push_back(r);
        deaths.push_back(d);
//        std::cout << c << ", " << r << ", " << d << std::endl;
      }
      in.close();
  }

  inline int getNumDays() const { return confirmed.size(); }
  inline int getNumActive(int day) const { return confirmed[day] - recovered[day] - deaths[day]; }

  int N = 0; //population
  std::vector<int> confirmed;
  std::vector<int> recovered;
  std::vector<int> deaths;
};

Vector getQuarantineInputVector(const std::string& filepath);

#endif
