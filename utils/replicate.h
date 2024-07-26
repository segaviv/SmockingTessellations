#pragma once

#include <Eigen/Eigen>
#include "vecops.h"
#include <vector>

namespace replicate {

// Replicate a vector while shifting each element by a certain amount on each
// replication.
// res[i] = v[i % v.size()] + shift * (i / v.size())
template<typename T, typename G>
std::vector<T> vec_with_shift(const std::vector<T>& v, int times, G shift) {
  std::vector<T> res(v.size() * times);
  for (int i = 0; i < times; ++i) {
    for (int j = 0; j < v.size(); ++j) {
      res[i * v.size() + j] = v[j] + shift * i;
    }
  }
  return res;
}

template<typename T, typename G>
std::vector<T> vec_with_shift(const std::vector<T>& v, int start, int end, G shift) {
  int times = end - start;
  std::vector<T> res(v.size() * times);
  for (int i = start; i < end; ++i) {
    for (int j = 0; j < v.size(); ++j) {
      res[(i - start) * v.size() + j] = v[j] + shift * i;
    }
  }
  return res;
}

}