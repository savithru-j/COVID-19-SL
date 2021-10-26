#ifndef OBSERVEDPOPULATION_H
#define OBSERVEDPOPULATION_H

#include <array>
#include <fstream>

#include "linearalgebra/Vector.h"
#include "utils/ErrorHandler.h"

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

  void smooth(const int T_smooth);

  int N = 0; //population
  std::vector<int> confirmed;
  std::vector<int> recovered;
  std::vector<int> deaths;
};

Vector<double> getQuarantineInputVector(const std::string& filepath);
Vector<double> getDailyVaccinations(const std::string& filepath, const int T_smooth);

#endif
