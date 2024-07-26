#include "seamless_arap.h"
#include "utils/Hmesh.h"
#include <chrono>
#include <igl/svd3x3.h>
#include <utils/conversions.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace visualization {
SeamlessARAP::SeamlessARAP(const utils::Hmesh &mesh,
                           const std::unordered_map<int, int> &stitch_verts)
    : mesh(mesh), stitch_verts(stitch_verts) {}

SeamlessARAP::SeamlessARAP(const SeamlessARAP &other)
    : one_rings(other.one_rings), old_to_new(other.old_to_new),
      new_to_old(other.new_to_old), mesh(other.mesh) {
  objective = std::make_unique<Optiz::QuadraticObjectiveD>(*other.objective);
}

SeamlessARAP &SeamlessARAP::set_objective() {
  build_one_rings();
  objective = std::make_unique<Optiz::QuadraticObjectiveD>(new_V.rows(), 3);
  // Add the ARAP energy.
  objective->add_weighted_equations(
      one_rings.size(), [&](int v, auto &x, const auto &add_eq) {
        auto &ring = one_rings[v];
        for (int ind = 0; ind < ring.v0.size(); ind++) {
          int i = ring.v0[ind], j = ring.v1[ind];
          double w_ij = ring.cot_weights[ind];
          // w_ij |(xj - xi) - rhs|^2.
          add_eq(w_ij, x(j) - x(i));
        }
      });
  return *this;
}

SeamlessARAP &
SeamlessARAP::set_fixed_indices(const std::vector<int> &fixed_indices) {
  objective->set_known_indices(
      Eigen::VectorXi::Map(fixed_indices.data(), fixed_indices.size()));
  return *this;
}

SeamlessARAP &
SeamlessARAP::set_fixed_indices(const Eigen::VectorXi &fixed_indices) {
  objective->set_known_indices(fixed_indices);
  return *this;
}

SeamlessARAP &SeamlessARAP::set_fixed(const std::vector<int> &indices,
                                      const Eigen::MatrixXd &positions) {
  std::set<int> used;
  std::vector<int> slice;
  std::vector<int> new_indices;
  slice.reserve(indices.size());
  new_indices.reserve(indices.size());
  for (int i = 0; i < indices.size(); i++) {
    if (!used.insert(old_to_new[indices[i]]).second)
      continue;
    slice.push_back(i);
    new_indices.push_back(old_to_new[indices[i]]);
  }
  objective->set_knowns(convert::eigen_vec_map(new_indices),
                        positions(slice, Eigen::all));
  return *this;
}

SeamlessARAP &SeamlessARAP::precompute() {
  objective->prefactor();
  return *this;
}

static Eigen::Matrix3d closest_rotation(const Eigen::Matrix3d &mat) {
  Eigen::Matrix3d U, V;
  Eigen::Vector3d S;
  igl::svd3x3(mat, U, S, V);
  // IGL's svd3x3 guarantees that U, V are rotation matrices.
  return U * V.transpose();
}

Eigen::MatrixXd SeamlessARAP::solve(const Eigen::MatrixXd &X, bool zero_rhs) {
  Eigen::MatrixXd rhs(objective->get_quadratic_term(0).A.rows(), 3);
  if (zero_rhs) {
    rhs.setZero();
    return objective->solve(rhs);
  }
// Calculate the RHS.
#pragma omp parallel for
  for (auto &one_ring : one_rings) {
    auto edges = X(one_ring.v1, Eigen::all) - X(one_ring.v0, Eigen::all);
    auto old_edges = one_ring.original_edges();
    // Find best fitting rotation from original to current edges.
    Eigen::Matrix3d R = closest_rotation(old_edges.transpose() * edges);
    // Add the rotated original edges to the RHS.
    one_ring.edges(rhs) = old_edges * R;
  }
  return objective->solve(rhs);
}

SeamlessARAP &
SeamlessARAP::set_fixed_positions(const Eigen::MatrixXd &fixed_positions) {
  objective->set_known_values(fixed_positions);
  return *this;
}

void SeamlessARAP::build_one_rings() {
  // Some vertices are merged, so we need to remap the vertices.
  build_verts_mapping();
  auto new_v = [&](int v) { return old_to_new[v]; };
  auto old_v = [&](int v) { return new_to_old[v]; };

  one_rings.resize(new_V.rows());
  auto add_edge = [&](int dest, int v0, int v1, double cot_weight) {
    int new_dest = new_v(dest);
    one_rings[new_dest].v0.push_back(new_v(v0));
    one_rings[new_dest].v1.push_back(new_v(v1));
    one_rings[new_dest].edges_vec.push_back(mesh.V.row(v1) - mesh.V.row(v0));
    // If we remove 'dest' and keep 'old_v(new_dest)', apply the seam rotation.
    if (old_v(new_dest) != dest) {
      one_rings[new_dest].edges_vec.back().head<2>() *=
          Eigen::Rotation2D(rotation_angles[dest])
              .toRotationMatrix()
              .transpose();
    }
    one_rings[new_dest].cot_weights.push_back(cot_weight);
  };

  for (auto &e : mesh.edges) {
    if (e.twin() && e.vi > e.next()->vi) {
      continue; // Don't add the same edge twice.
    }
    double cot_weight = e.cot_theta();
    double twin_cot_weight = (e.twin() ? e.twin()->cot_theta() : 0);
    int v0 = e.vi, v1 = e.next()->vi;
    // Add the edge to the neighborhoods of the triangle vertices.
    add_edge(e.vi, v0, v1, cot_weight + twin_cot_weight);
    add_edge(e.next()->vi, v0, v1, cot_weight + twin_cot_weight);
    add_edge(e.prev()->vi, v0, v1, cot_weight);
  }
  int current_row = 0;
  for (int v = 0; v < one_rings.size(); v++) {
    one_rings[v].start_row = current_row;
    current_row += one_rings[v].num_edges();
  }
}

double SeamlessARAP::find_seam_rotation(int v0, int v1) {
  // Helper funcs.
  auto are_dups = [&](int v0, int v1) {
    auto it = stitch_verts.find(v0);
    return it != stitch_verts.end() && it->second == v1;
  };
  auto angle_between = [&](int p0, int p1, int q0, int q1) {
    Eigen::VectorXd e0 = mesh.V.row(p1) - mesh.V.row(p0),
                    e1 = mesh.V.row(q1) - mesh.V.row(q0);
    return std::atan2(e1.y(), e1.x()) - std::atan2(e0.y(), e0.x());
  };
  // Find two stitching edges and calculate rotation.
  int next_v0 = mesh.next_boundary_vert(v0),
      prev_v1 = mesh.prev_boundary_vert(v1);
  if (are_dups(next_v0, prev_v1))
    return angle_between(v0, next_v0, v1, prev_v1);

  int prev_v0 = mesh.prev_boundary_vert(v0),
      next_v1 = mesh.next_boundary_vert(v1);
  if (are_dups(prev_v0, next_v1))
    return angle_between(v0, prev_v0, v1, next_v1);

  std::cout << "Error: Failed to find seam rotation." << std::endl;
  return 0;
}

void SeamlessARAP::build_verts_mapping() {
  old_to_new = std::vector<int>(mesh.nv(), -1);
  new_to_old.clear();
  rotation_angles = std::vector<double>(mesh.nv(), 0);

  for (int i = 0; i < mesh.nv(); i++) {
    if (old_to_new[i] != -1)
      continue;
    old_to_new[i] = new_to_old.size();
    new_to_old.push_back(i);
    auto it = stitch_verts.find(i);
    if (it != stitch_verts.end()) {
      old_to_new[it->second] = old_to_new[i];
      rotation_angles[it->second] = find_seam_rotation(it->second, i);
    }
  }
  // Update new_V and new_F.
  new_V = mesh.V(new_to_old, Eigen::all);
  new_F = mesh.F;
  for (auto &face : new_F)
    for (auto &v : face)
      v = old_to_new[v];
}

} // namespace visualization