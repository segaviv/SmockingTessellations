#pragma once

#include <Eigen/Eigen>
#include <unordered_map>
#ifndef EMSCRIPTEN
#include <glm/glm.hpp>
#endif
#include <vector>

namespace convert {

template <typename T>
std::vector<Eigen::VectorX<T>> to_vec_mat(const Eigen::MatrixX<T> &m) {
  std::vector<Eigen::VectorX<T>> v(m.rows());
  for (int i = 0; i < m.rows(); ++i) {
    v[i] = m.row(i);
  }
  return v;
}

template <typename T, int R, int C>
Eigen::MatrixX<T> to_eig_mat(const std::vector<Eigen::Matrix<T, R, C>> &v) {
  Eigen::MatrixX<T> m(v.size(), v[0].size());
  for (int i = 0; i < v.size(); ++i) {
    m.row(i) = v[i];
  }
  return m;
}
template <typename T>
Eigen::MatrixX<T> to_eig_mat(const std::vector<std::vector<T>> &mat) {
  Eigen::MatrixX<T> m(mat.size(), mat[0].size());
  for (int i = 0; i < mat.size(); ++i) {
    for (int j = 0; j < mat[0].size(); ++j) {
      m(i, j) = mat[i][j];
    }
  }
  return m;
}

template <typename T> Eigen::MatrixX<T> to_3d(const Eigen::MatrixX<T> &m) {
  if (m.cols() == 3) {
    return m;
  }
  Eigen::MatrixX<T> m3d(m.rows(), 3);
  m3d << m, Eigen::MatrixX<T>::Zero(m.rows(), 1);
  return m3d;
}

template <typename T, typename Container> T to(const Container &v) {
  return T(v.begin(), v.end());
}

template <typename T, typename G>
std::vector<T> to_vector(const std::unordered_map<T, G> &m) {
  std::vector<T> v;
  v.reserve(m.size());
  for (const auto &p : m)
    v.push_back(p.first);
  return v;
}

#ifndef EMSCRIPTEN
inline glm::vec3 to_glm(const Eigen::Vector3d &v) {
  return glm::vec3(v.x(), v.y(), v.z());
}
#endif

template<typename T>
auto eigen_vec_map(const std::vector<T>& vec) {
  return Eigen::VectorX<T>::Map(vec.data(), vec.size());
}

} // namespace convert