#pragma  once

#include <vector>
#include <set>

namespace slice {

  template<typename T, typename IntContainer>
  std::vector<T> vec(const std::vector<T>& v, const IntContainer& indices) {
    std::vector<T> sliced(indices.size());
    int cur = 0;
    for (int i : indices) {
      sliced[cur++] = v[i];
    }
    return sliced;
  }

  template<typename Func>
  std::vector<int> indices(int n, const Func& func) {
    std::vector<int> res;
    for (int i = 0; i < n; i++) {
      if (func(i))
        res.push_back(i);
    }
    return res;
  }
}