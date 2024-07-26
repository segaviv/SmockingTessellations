#pragma once

#include "utils/Hmesh.h"
#include <Eigen/Eigen>
#include <Optiz/Linear/QuadraticObjectiveD.h>
#include <memory>
#include <unordered_map>
#include <vector>
namespace visualization {

class SeamlessARAP {
public:
  SeamlessARAP(const utils::Hmesh &mesh,
               const std::unordered_map<int, int> &stitch_verts =
                   std::unordered_map<int, int>());

  SeamlessARAP(const SeamlessARAP &other);

  SeamlessARAP &set_objective();
  SeamlessARAP &set_fixed_indices(const std::vector<int> &fixed_indices);
  SeamlessARAP &set_fixed_indices(const Eigen::VectorXi &fixed_indices);
  SeamlessARAP &precompute();
  Eigen::MatrixXd solve(const Eigen::MatrixXd &X, bool zero_rhs = false);
  SeamlessARAP &set_fixed_positions(const Eigen::MatrixXd &fixed_positions);

  SeamlessARAP &set_fixed(const std::vector<int> &indices,
                          const Eigen::MatrixXd &positions);

  struct VertexOneRing {
    VertexOneRing() {
      v0.reserve(16);
      v1.reserve(16);
      edges_vec.reserve(16);
      cot_weights.reserve(16);
    }
    // The spokes and rim edge indices.
    std::vector<int> v0, v1;
    // Corresponding cotangent weights.
    std::vector<double> cot_weights;
    std::vector<Eigen::RowVector3d> edges_vec;
    // The original edges (V(v0, :) - V(v1, :)).
    auto original_edges() const {
      return Eigen::Matrix<double, -1, 3, Eigen::RowMajor>::Map(
          (double *)edges_vec.data(), v0.size(), 3);
    };
    // [start_row, start_row + num_edges()) contain the linear equations
    // corresponding to this one ring.
    int start_row;
    int num_edges() { return v0.size(); }
    Eigen::Block<Eigen::Matrix<double, -1, -1>> edges(Eigen::MatrixXd& rhs) {
      return rhs.middleRows(start_row, num_edges());
    }
  };

  Eigen::MatrixXd verts() const { return new_V; }
  std::vector<std::vector<int>> faces() const { return new_F; }

private:
  void build_one_rings();
  void build_verts_mapping();
  double find_seam_rotation(int v0, int v1);

  std::vector<VertexOneRing> one_rings;
  utils::Hmesh mesh;
  // Stitching information.
  std::unordered_map<int, int> stitch_verts;
  std::vector<double> rotation_angles;
  std::vector<int> old_to_new, new_to_old;

  // If we stitch vertices, some of them will be merged.
  Eigen::MatrixXd new_V;
  std::vector<std::vector<int>> new_F;

  std::unique_ptr<Optiz::QuadraticObjectiveD> objective;
};
} // namespace visualization