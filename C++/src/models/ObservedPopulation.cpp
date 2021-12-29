#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "ObservedPopulation.h"

void
ObservedPopulation::smooth(const int T_smooth)
{
  assert(T_smooth % 2 == 1);
  const int r = (T_smooth-1)/2;

  auto confirmed_tmp = confirmed;
  auto recovered_tmp = recovered;
  auto deaths_tmp = deaths;

  for (int i = 0; i < (int) confirmed.size(); ++i)
  {
    int js = std::max(i-r, 0);
    int je = std::min(i+r, (int)confirmed.size()-1);
    int r_i = std::min(i-js, je-i);
    double sum_c = 0, sum_r = 0, sum_d = 0;
    double denom = 2*r_i + 1;
    for (int j = i-r_i; j <= i+r_i; ++j)
    {
      sum_c += confirmed_tmp[j];
      sum_r += recovered_tmp[j];
      sum_d += deaths_tmp[j];
    }
    confirmed[i] = std::round(sum_c / denom);
    recovered[i] = std::round(sum_r / denom);
    deaths[i]    = std::round(sum_d / denom);
  }
}

Vector<double>
getQuarantineInputVector(const std::string& filepath)
{
  std::ifstream in(filepath);
  if (!in)
    return Vector<double>(); //return empty vector

  std::vector<int> indices;
  std::vector<double> values;

  int ind, max_ind = -1;
  while (in >> ind)
  {
    double val;
    in >> val;
    indices.push_back(ind);
    values.push_back(val);
    max_ind = std::max(max_ind, ind);
  }
  in.close();

  Vector<double> vec(max_ind+1, 0.0);
  for (std::size_t i = 0; i < indices.size(); ++i)
    vec[indices[i]] = values[i];

  return vec;
}

Vector<double>
readDataVectorFromFile(const std::string& filepath, const int T_smooth)
{
  std::ifstream in(filepath);
  if (!in)
    return Vector<double>(); //return empty vector

  Vector<double> original_data;
  double val;
  while (in >> val)
    original_data.push_back(val);
  in.close();

  if (T_smooth == 0)
    return original_data;

  assert(T_smooth % 2 == 1);
  const int r = (T_smooth-1)/2;

  Vector<double> smoothed_data(original_data.size(), 0.0);
  for (int i = 0; i < (int) original_data.size(); ++i)
  {
    int js = std::max(i-r, 0);
    int je = std::min(i+r, (int)original_data.size()-1);
    int r_i = std::min(i-js, je-i);
    double sum_v = 0;
    double denom = 2*r_i + 1;
    for (int j = i-r_i; j <= i+r_i; ++j)
      sum_v += original_data[j];
    smoothed_data[i] = std::round(sum_v / denom);
  }

  //Compute diff to get daily values
//  for (int i = (int)original_data.size()-1; i > 0; --i)
//    daily_vaccinated_smooth[i] -= daily_vaccinated_smooth[i-1];
//  daily_vaccinated_smooth[0] = 0;

  return smoothed_data;
}
