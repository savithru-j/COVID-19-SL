#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
#include <array>
#include "ErrorHandler.h"

template<class T = double>
class Vector
{
public:

  Vector() = default;
  Vector(int m) : data_(m) {}
  Vector(int m, const T& val) : data_(m, val) {}

  Vector(const std::vector<T>& v) : data_(v) {}

  inline int m() const { return data_.size(); }
  inline std::size_t size() const { return data_.size(); }

  inline const T& operator()(int i) const { return data_[i]; }
  inline T& operator()(int i) { return data_[i]; }

  inline const T& operator[](int i) const { return data_[i]; }
  inline T& operator[](int i) { return data_[i]; }

  inline const T& back() const { return data_.back(); }

  inline typename std::vector<T>::iterator begin() { return data_.begin(); }
  inline typename std::vector<T>::iterator end() { return data_.end(); }

  inline void resize(int m, const T& val = 0.0) { data_.resize(m, val); }

  inline void clear() { data_.clear(); }

  inline void push_back(const T& val) { data_.push_back(val); }

  inline void insert(typename std::vector<T>::const_iterator pos, const T& val)
  {
    data_.insert(pos, val);
  }

  inline const std::vector<T>& getDataVector() const { return data_; }
  inline std::vector<T>& getDataVector() { return data_; }

protected:
  std::vector<T> data_;
};

template<class T = double>
class Matrix
{
public:

  Matrix() = default;
  Matrix(int m, int n, double val = 0.0)
    : m_(m), n_(n), data_(m*n, val) {}

  inline int m() const { return m_; }
  inline int n() const { return n_; }

  inline const T& operator()(int i, int j) const { return data_[n_*i + j]; }
  inline T& operator()(int i, int j) { return data_[n_*i + j]; }

protected:
  int m_ = 0, n_ = 0;
  std::vector<T> data_;
};

template<class T>
inline Vector<T> operator*(const Matrix<T>& A, const Vector<T>& x)
{
  if (A.n() != x.m())
    throwError("Dimension mismatch for matrix vector multiplication!");

  int m = A.m();
  int n = A.n();
  Vector<T> y(m, 0.0);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      y(i) += A(i,j)*x(j);

  return y;
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
{
  for (int i = 0; i < v.m()-1; ++i)
    os << v[i] << ", ";
  os << v.back();
  return os;
}

template<typename T, std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& v)
{
  for (std::size_t i = 0; i < N-1; ++i)
    os << v[i] << ", ";
  os << v.back();
  return os;
}

#endif
