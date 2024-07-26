#include "UnitSmockingPattern.h"

namespace smocking {

UnitSmockingPattern UnitSmockingPattern::ARROW = UnitSmockingPattern(
    5, 3, {{{0, 2}, {1, 1}, {2, 2}}, {{2, 1}, {3, 0}, {4, 1}}});

UnitSmockingPattern UnitSmockingPattern::LEAF = UnitSmockingPattern(
    7, 3,
    {{{0, 1}, {1, 0}}, {{1, 1}, {2, 2}}, {{3, 0}, {4, 1}}, {{4, 2}, {5, 1}}});

UnitSmockingPattern UnitSmockingPattern::BRAID = UnitSmockingPattern(
    4, 3,
    {{{0, 1}, {1, 2}}, {{1, 1}, {2, 0}}});

    UnitSmockingPattern UnitSmockingPattern::HEART = UnitSmockingPattern(
    7, 3,
    {{{0, 1}, {1, 0}}, 
    {{1, 1}, {2, 2}},
    {{3, 2}, {4, 1}},
    {{4, 0}, {5, 1}}
    });

UnitSmockingPattern UnitSmockingPattern::BOX = UnitSmockingPattern(
    5, 7,
    {{{1, 1}, {0, 2}}, {{1, 2}, {0, 3}},
    {{2, 1}, {3, 2}}, {{2, 2}, {3, 3}},
    {{0, 4}, {1, 5}}, {{0, 5}, {1, 6}},
    {{3, 4}, {2, 5}}, {{3, 5}, {2, 6}}});
    UnitSmockingPattern UnitSmockingPattern::BRICK = UnitSmockingPattern(
    5, 5,
    {{{0, 0}, {1, 1}, {0, 2}}, 
    {{1, 2}, {0, 3}, {1,4}},
    {{3, 0}, {2, 1}, {3, 2}},
    {{2, 2}, {3, 3}, {2,4}}});

UnitSmockingPattern UnitSmockingPattern::TWISTED_SQUARE = UnitSmockingPattern(
    4, 4,
    {{{0, 1}, {1, 1}}, {{2, 1}, {2, 2}},
    {{0, 2}, {0, 3}}, {{1, 3}, {2, 3}}});


UnitSmockingPattern::UnitSmockingPattern(int nx, int ny) : _nx(nx), _ny(ny) {
  _verts.resize(nx * ny, 2);
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      _verts.row(coords_to_index(x, y)) = Eigen::Vector2d(x, y);
    }
  }
}

UnitSmockingPattern::UnitSmockingPattern(
    int nx, int ny,
    const std::vector<std::vector<Eigen::Vector2i>> &stitching_lines)
    : UnitSmockingPattern(nx, ny) {
  _stitching_lines = stitching_lines;
}

void UnitSmockingPattern::add_stitching_line(
    const std::vector<Eigen::Vector2i> &line) {
  _stitching_lines.push_back(line);
}

std::vector<Eigen::Vector2i> UnitSmockingPattern::get_stitching_pairs() const {
  std::vector<Eigen::Vector2i> pairs;
  for (const auto &line : _stitching_lines) {
    for (int i = 0; i < line.size() - 1; i++) {
      pairs.push_back({coords_to_index(line[i]), coords_to_index(line[i + 1])});
    }
  }
  return pairs;
}

Eigen::MatrixXd UnitSmockingPattern::get_coord_to_sl_index() const {
  Eigen::MatrixXd coord_to_sl_index = Eigen::MatrixXd::Constant(_nx, _ny, -1);
  for (int i = 0; i < _stitching_lines.size(); i++) {
    for (const auto &coords : _stitching_lines[i]) {
      coord_to_sl_index(coords.x(), coords.y()) = i;
    }
  }
  return coord_to_sl_index;
}

} // namespace smocking