#pragma once

#include "UnitSmockingPattern.h"
#include <utils/Hmesh.h>
#include <vector>

namespace smocking {
class Tangram {
public:
  enum FaceClassification {
    NOT_SET = 0,
    UNDERLAY = 1,
    PLEAT = 2,
  };

  Tangram(const UnitSmockingPattern &pattern);

  std::vector<Eigen::Vector2i> get_underlay_edges() const {
    return _underlay_edges;
  }
  const std::vector<std::vector<int>>& get_triangles() const {
    return _triangulated_faces;
  }
  const std::vector<FaceClassification>& get_face_classification() const {
    return _face_classification;
  }

  const UnitSmockingPattern &get_pattern() const { return _pattern; }

  bool is_underlay_vertex(int index) const;
  // Stitching/underlay edge.
  bool is_underlay_edge(const Eigen::Vector2i &coord1,
                        const Eigen::Vector2i &coord2) const;
  bool is_underlay_edge(int index1, int index2) const;
  bool is_stitching_edge(const Eigen::Vector2i &coord1,
                         const Eigen::Vector2i &coord2) const;
  bool is_stitching_edge(int index1, int index2) const;
  int index_to_sl(int index) const;

  // Get the optimized vertices positions. We hire Gary to control
  // the stitching scale (Gary = 1.0 - no stitching, Gary = 0.0 - fully
  // stitched).
  Eigen::MatrixXd compute_vertices_positions(double Gary);

  // Set the stitching scale and update the vertices positions.
  void set_stitching_scale(double Gary);

  // Returns the basis for the periodicity of the pattern.
  // B = [v1; v2], P_0 - unit pattern. The tiling is given by
  // P = P_0 + v1 * x + v2 * y.
  Eigen::Matrix2d compute_periodic_shift(const Eigen::MatrixXd& verts_positions) const;
  Eigen::Matrix2d compute_periodic_basis() const;

  // Returns a subset of the mesh representing the tangram face that
  // contains the start_triangle.
  utils::Hmesh get_tangram_face(int start_triangle);

  // Mean value coordinates of a pleat vertex.
  struct PleatMVC {
    struct SurroundingVertex {
      // Index into the vertices.
      int index;
      // Weight.
      double weight;
      Eigen::RowVector2d shift;
    };
    int pleat_index;
    std::vector<SurroundingVertex> coords;
    void normalize_coords() {
      double sum = 0.0;
      for (auto &coord : coords) {
        sum += coord.weight;
      }
      for (auto &coord : coords) {
        coord.weight /= sum;
      }
    }
  };
  std::vector<PleatMVC> get_pleat_mvc() const { return _pleat_mvc; }

private:
double stitching_lines_distance(int sl1, int sl2) const;

public:
  void apply_periodic_bc(Eigen::MatrixXd &coord_to_sl);
  void apply_periodic_bc(utils::Hmesh &mesh);
  void init_underlay_edges();
  void triangulate();
  void classify_triangles();
  void compute_pleats_mvc();

  // Triangulated grid with face classification.
  std::vector<std::vector<int>> _triangulated_faces;
  std::vector<FaceClassification> _face_classification;
  utils::Hmesh tangram_mesh;

  std::vector<Eigen::Vector2i> _underlay_edges;
  UnitSmockingPattern _pattern;

  std::vector<PleatMVC> _pleat_mvc;
  // Indices of the vertices that are pleats or underlays.
  std::vector<int> _pleat_indices;
  std::vector<int> _underlay_indices;
  std::vector<int> _index_to_pleat_mvc_index;

  // Sewing line index for each coordinate in the grid (-1 if no sewing line).
  Eigen::MatrixXd _coord_to_sl;
  std::vector<std::vector<Eigen::Vector2i>> _stitching_lines;

  // The stitching scale and the vertices positions computed by optimization.
  double _stitching_scale = 1.0;
  Eigen::MatrixXd _vertices_positions;
};
} // namespace smocking
