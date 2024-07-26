#pragma once
#include <Eigen/Eigen>
#include <vector>

template<typename T>
std::vector<T> operator+(const std::vector<T>& v, T shift) {
  std::vector<T> res(v.size());
  for (int i = 0; i < v.size(); ++i) {
    res[i] = v[i] + shift;
  }
  return res;
}

template<typename T, int Rows, int Cols>
Eigen::Matrix<T, Rows, Cols> operator+(const Eigen::Matrix<T, Rows, Cols>& m, T shift) {
  return m.array() + shift;
}