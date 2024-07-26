#pragma once

#include <Eigen/Core>

namespace utils {
  std::pair<int, Eigen::Vector3d> point_in_param(const Eigen::MatrixXd &UV, const std::vector<std::vector<int>> &F, const Eigen::Vector2d &p);
}