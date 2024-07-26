#pragma once

#include <Eigen/Eigen>
#include <functional>
#include <vector>

namespace construct {

template<typename T>
std::vector<T> vec(int num_elements, std::function<T(int)> elem_func) {
  std::vector<T> v(num_elements);
  for (int i = 0; i < num_elements; ++i) {
    v[i] = elem_func(i);
  }
  return v;
}

template<typename T>
std::vector<T> linspace(T start, T end, T step = T(1)) {
  std::vector<T> v;
  for (T i = start; i < end; i += step) {
    v.push_back(i);
  }
  return v;
}

template<typename T>
std::vector<T> linspace(T end) {
  return linspace(0, end);
}


}