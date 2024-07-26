#pragma once

#include <directional/CartesianField.h>
#include <smocking_pattern/Tangram.h>
#include <unordered_map>
#include <vector>

namespace param {

using FaceClassification = smocking::Tangram::FaceClassification;

struct PullbackResult {
  // A pulled back vertex.
  struct Vertex {
    int unit_pattern_index;    // Index of the vertex in the unit pattern.
    Eigen::RowVector2d repeat; // Index of the vertex in the repeated pattern.

    Eigen::Vector2d open_pos; // The position in the open configuration.

    int uv_face;
    Eigen::Vector3d bary_coords;
    bool is_underlay;
    bool is_inside;
    bool is_boundary_pleat;
  };
  std::vector<Vertex> vertices;
  std::vector<std::vector<int>> faces;

  // The indices of the vertices that are inside the tangram.
  std::vector<int> inside_verts_indices;
  std::vector<int> inside_underlay_indices;
  std::vector<std::vector<int>> reduced_faces;
  // Convert from the original index to the reduced index.
  std::vector<int> orig_to_inside;

  std::vector<FaceClassification> face_classification;
  std::vector<Eigen::Vector2i> underlay_edges;
  std::vector<Eigen::Vector2i> stitching_edges;
  // Mean value coordinates of the pleat vertices.
  std::vector<smocking::Tangram::PleatMVC> pleat_mvcs;

  // Stitching vertices (first should be stitched to all seconds).
  std::unordered_map<int, std::set<int>> duplicate_verts;

  // Hacky, but whatever...
  std::function<Eigen::MatrixXd(double)> get_open_configuration;
  std::function<Eigen::MatrixXd()> get_closed_configuration;

  // Get the open (flat) configuration.
  Eigen::MatrixXd get_open_configuration_no_symmetry(double opening) const;
  // Get the closed (3d) configuration.
  Eigen::MatrixXd get_closed_configuration_no_symmetry() const;
};

void pullback_no_rosy(double rotate_uv_angle = 0, double scale = 1.0,
                      const Eigen::Vector2d &shift = Eigen::Vector2d::Zero());

} // namespace param