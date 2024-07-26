#include "point_in_param.h"
#include "utils/Hmesh.h"

namespace utils {

  Eigen::Vector3d barycentric_coordinates(const Eigen::Vector2d &v0, const Eigen::Vector2d &v1, const Eigen::Vector2d &v2, const Eigen::Vector2d &p) {
    Eigen::Vector2d v0v1 = v1 - v0, v0v2 = v2 - v0, v0p = p - v0;
    double d00 = v0v1.dot(v0v1);
    double d01 = v0v1.dot(v0v2);
    double d11 = v0v2.dot(v0v2);
    double d20 = v0p.dot(v0v1);
    double d21 = v0p.dot(v0v2);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    return Eigen::Vector3d(1 - v - w, v, w);
  }

  std::pair<int, Eigen::Vector3d> point_in_param(const Eigen::MatrixXd &UV, const std::vector<std::vector<int>> &F, const Eigen::Vector2d &p) {
    // Find the triangle that contains the point p
    for (int i = 0; i < F.size(); i++) {
      Eigen::Vector2d v0 = UV.row(F[i][0]);
      Eigen::Vector2d v1 = UV.row(F[i][1]);
      Eigen::Vector2d v2 = UV.row(F[i][2]);
      Eigen::Vector3d bary = utils::barycentric_coordinates(v0, v1, v2, p);
      if (bary(0) >= 0 && bary(1) >= 0 && bary(2) >= 0) {
        return std::make_pair(i, bary);
      }
    }
    return std::make_pair(-1, Eigen::Vector3d::Zero());
  }
}