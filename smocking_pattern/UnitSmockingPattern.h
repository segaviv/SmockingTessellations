#pragma once

#include <Eigen/Eigen>
#include <vector>

namespace smocking {

class UnitSmockingPattern {
public:
  UnitSmockingPattern(int nx, int ny);
  UnitSmockingPattern(int nx, int ny, const std::vector<std::vector<Eigen::Vector2i>>& stitching_lines);

  void add_stitching_line(const std::vector<Eigen::Vector2i> &line);

  // Helper functions for converting between 2D coordinates and 1D indices.
  int coords_to_index(const Eigen::Vector2i &coords) const {
    return coords_to_index(coords.x(), coords.y());
  }
  int coords_to_index(int x, int y) const { return y * _nx + x; }
  Eigen::Vector2i index_to_coords(int index) const {
    return Eigen::Vector2i(index % _nx, index / _nx);
  }

  // Get pairs of vertices that should be stitched together.
  std::vector<Eigen::Vector2i> get_stitching_pairs() const;
  // Get a matrix that maps 2D coordinates to stitching line index.
  Eigen::MatrixXd get_coord_to_sl_index() const;

  int nx() const { return _nx; }
  int ny() const { return _ny; }
  std::vector<std::vector<Eigen::Vector2i>> stitching_lines() const { return _stitching_lines; }
  int num_sl() const { return _stitching_lines.size(); }
  const Eigen::MatrixXd& verts() const { return _verts; }

  static UnitSmockingPattern ARROW;
  static UnitSmockingPattern LEAF;
  static UnitSmockingPattern BRAID;
  static UnitSmockingPattern HEART;
  static UnitSmockingPattern BOX;
  static UnitSmockingPattern BRICK;
  static UnitSmockingPattern TWISTED_SQUARE;

private:
  int _nx;
  int _ny;
  std::vector<std::vector<Eigen::Vector2i>> _stitching_lines;
  Eigen::MatrixXd _verts;
};

} // namespace smocking
