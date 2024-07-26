#include "Tangram.h"
#include "Optiz/NewtonSolver/Problem.h"
#include "smocking_pattern/UnitSmockingPattern.h"
#include "utils/Hmesh.h"
#include "utils/conversions.h"
#include "utils/slicing.h"
#include <Optiz/Problem.h>
#include <cassert>
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

namespace smocking {

Tangram::Tangram(const UnitSmockingPattern &pattern) : _pattern(pattern) {
  _vertices_positions = _pattern.verts();
  _coord_to_sl = pattern.get_coord_to_sl_index();
  _stitching_lines = pattern.stitching_lines();
  apply_periodic_bc(_coord_to_sl);
  for (int i = 0; i < _coord_to_sl.size(); i++) {
    if (_coord_to_sl(i) == -1) {
      _pleat_indices.push_back(i);
    } else {
      _underlay_indices.push_back(i);
    }
  }
  // Find underlay edges.
  init_underlay_edges();
  // Triangulate the grid faces according to the underlay edges.
  triangulate();
  // Compute the MVC for pleats.
  compute_pleats_mvc();
}

void Tangram::apply_periodic_bc(Eigen::MatrixXd &coord_to_sl) {
  int nx = _pattern.nx(), ny = _pattern.ny();
  _stitching_lines.resize(7 * _pattern.num_sl());
  for (int y = 0; y < ny; y++) {
    if (coord_to_sl(nx - 1, y) == -1 && coord_to_sl(0, y) != -1) {
      int new_sl = coord_to_sl(0, y) + _pattern.num_sl();
      coord_to_sl(nx - 1, y) = new_sl;
      _stitching_lines[new_sl].push_back({nx - 1, y});
    }
    if (coord_to_sl(0, y) == -1 && coord_to_sl(nx - 1, y) != -1) {
      int new_sl = coord_to_sl(nx - 1, y) + 2 * _pattern.num_sl();
      coord_to_sl(0, y) = new_sl;
      _stitching_lines[new_sl].push_back({0, y});
    }
  }
  for (int x = 0; x < nx; x++) {
    if (coord_to_sl(x, ny - 1) == -1 && coord_to_sl(x, 0) != -1) {
      int new_sl = coord_to_sl(x, 0) + 3 * _pattern.num_sl();
      coord_to_sl(x, ny - 1) = new_sl;
      _stitching_lines[new_sl].push_back({x, ny - 1});
    }
    if (coord_to_sl(x, 0) == -1 && coord_to_sl(x, ny - 1) != -1) {
      int new_sl = coord_to_sl(x, ny - 1) + 4 * _pattern.num_sl();
      coord_to_sl(x, 0) = new_sl;
      _stitching_lines[new_sl].push_back({x, 0});
    }
  }
}

void Tangram::apply_periodic_bc(utils::Hmesh &mesh) {
  for (int i = 0; i < mesh.edges.size(); i++) {
    if (!mesh.edges[i].is_boundary())
      continue;
    int v0 = mesh.edges[i].vi, v1 = mesh.edges[i].next()->vi;
    Eigen::Vector2i coord0 = _pattern.index_to_coords(v0),
                    coord1 = _pattern.index_to_coords(v1);
    if (coord0.x() == 0 && coord1.x() == 0) {
      coord0.x() = _pattern.nx() - 1;
      coord1.x() = _pattern.nx() - 1;
    } else if (coord0.x() == _pattern.nx() - 1 &&
               coord1.x() == _pattern.nx() - 1) {
      coord0.x() = 0;
      coord1.x() = 0;
    } else if (coord0.y() == 0 && coord1.y() == 0) {
      coord0.y() = _pattern.ny() - 1;
      coord1.y() = _pattern.ny() - 1;
    } else if (coord0.y() == _pattern.ny() - 1 &&
               coord1.y() == _pattern.ny() - 1) {
      coord0.y() = 0;
      coord1.y() = 0;
    } else {
      continue;
    }
    // Other edge vertices.
    int v0_new = _pattern.coords_to_index(coord0),
        v1_new = _pattern.coords_to_index(coord1);
    if (!mesh.verts_edges[v1_new].count(v0_new))
      continue;
    auto edge = mesh.verts_edges[v1_new][v0_new];
    mesh.edges[i].ti = edge;
    mesh.edges[edge].ti = i;
  }
}

void Tangram::init_underlay_edges() {
  // Check for adjacent underlay vertices from different stitching lines.
  static std::vector<Eigen::Vector2i> dirs = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
  int nx = _pattern.nx(), ny = _pattern.ny();
  for (int i = 0; i < _coord_to_sl.size(); i++) {
    Eigen::Vector2i coord = _pattern.index_to_coords(i);
    for (auto &dir : dirs) {
      Eigen::Vector2i neighbor = coord + dir;
      if (neighbor.x() < 0 || neighbor.x() >= nx || neighbor.y() < 0 ||
          neighbor.y() >= ny)
        continue;

      // int sl1 = _coord_to_sl(coord.x(), coord.y());
      // int sl2 = _coord_to_sl(neighbor.x(), neighbor.y());
      // if (sl1 != sl2 && sl1 != -1 && sl2 != -1)
      // If both are underlay from different stitching lines, this is an
      // underlay edge.
      if (is_underlay_edge(coord, neighbor)) {
        _underlay_edges.push_back({_pattern.coords_to_index(coord),
                                   _pattern.coords_to_index(neighbor)});
      }
    }
  }
}

bool Tangram::is_underlay_vertex(int index) const {
  Eigen::Vector2i coord = _pattern.index_to_coords(index);
  return _coord_to_sl(coord.x(), coord.y()) != -1;
}

double Tangram::stitching_lines_distance(int sl1, int sl2) const {
  double dist = std::numeric_limits<double>::max();
  std::cout << "sl1: " << sl1 << " sl2: " << sl2 << ", total: " << _stitching_lines.size() << "\n";
  for (const auto &sl : _stitching_lines[sl1]) {
    for (const auto &other_sl : _stitching_lines[sl2]) {
      dist = std::min(dist, (sl - other_sl).cast<double>().norm());
    }
  }
  return dist;
}

bool Tangram::is_underlay_edge(const Eigen::Vector2i &coord1,
                               const Eigen::Vector2i &coord2) const {
  int sl1 = _coord_to_sl(coord1.x(), coord1.y()),
      sl2 = _coord_to_sl(coord2.x(), coord2.y());
  return sl1 != -1 && sl2 != -1 && sl1 != sl2
  && abs((stitching_lines_distance(sl1, sl2) - (coord1-coord2).cast<double>().norm())) < 1e-6;
}

bool Tangram::is_underlay_edge(int index1, int index2) const {
  if (index1 == 0 && index2 == 2) {
    std::cout << "Checking edge 0-2\n";
    std::cout << "Coord1: " << _pattern.index_to_coords(index1).transpose()
              << "\n";
    std::cout << "Coord2: " << _pattern.index_to_coords(index2).transpose()
              << "\n";
  }
  return is_underlay_edge(_pattern.index_to_coords(index1),
                          _pattern.index_to_coords(index2));
}

bool Tangram::is_stitching_edge(const Eigen::Vector2i &coord1,
                                const Eigen::Vector2i &coord2) const {
  int sl1 = _coord_to_sl(coord1.x(), coord1.y()),
      sl2 = _coord_to_sl(coord2.x(), coord2.y());
  return sl1 != -1 && sl2 != -1 && sl1 == sl2;
}

bool Tangram::is_stitching_edge(int index1, int index2) const {
  return is_stitching_edge(_pattern.index_to_coords(index1),
                           _pattern.index_to_coords(index2));
}

int Tangram::index_to_sl(int index) const {
  Eigen::Vector2i coord = _pattern.index_to_coords(index);
  return _coord_to_sl(coord.x(), coord.y());
}

void Tangram::triangulate() {
  // Helper functions.
  auto c2i = [&](const auto &coord) { return _pattern.coords_to_index(coord); };

  int nx = _pattern.nx(), ny = _pattern.ny();
  for (int x = 0; x < nx - 1; x++) {
    for (int y = 0; y < ny - 1; y++) {
      Eigen::Vector2i bl(x, y), br(x + 1, y), tr(x + 1, y + 1), tl(x, y + 1);
      // Determine how to triangulate the square (which diagonal to keep).
      if (is_underlay_edge(bl, tr) || is_stitching_edge(br, tl)) {
        _triangulated_faces.push_back({c2i(bl), c2i(br), c2i(tr)});
        _triangulated_faces.push_back({c2i(bl), c2i(tr), c2i(tl)});
      } else {
        _triangulated_faces.push_back({c2i(bl), c2i(br), c2i(tl)});
        _triangulated_faces.push_back({c2i(br), c2i(tr), c2i(tl)});
      }
    }
  }
  tangram_mesh = utils::Hmesh(_pattern.verts(), _triangulated_faces);
  apply_periodic_bc(tangram_mesh);
  classify_triangles();
}

// DFS from the given face until finding one with an underlay edge.
static std::pair<utils::Hmesh::Edge *, Eigen::RowVector2d>
get_bounding_underlay_edge(int start_triangle, Tangram &tangram) {
  std::set<int> visited;
  auto &mesh = tangram.tangram_mesh;
  std::vector<std::pair<int, Eigen::RowVector2d>> stack = {
      {start_triangle, {0, 0}}};
  while (!stack.empty()) {
    auto [face, shift] = stack.back();
    stack.pop_back();
    if (!visited.insert(face).second)
      continue;
    auto e = mesh.face(face).edge();
    for (int i = 0; i < 3; i++, e = e->next()) {
      if (tangram.is_underlay_edge(e->vi, e->next()->vi)) {
        return {e, shift};
      }
      // Account for periodic boundary conditions.
      auto new_shift =
          shift - mesh.V.row(e->twin()->vi) + mesh.V.row(e->next()->vi);
      stack.push_back({e->twin()->fi, new_shift});
    }
  }
  return {nullptr, {0, 0}};
}

// Get the next underlay edge bounding the underlay face.
static std::pair<utils::Hmesh::Edge *, Eigen::RowVector2d>
get_next_underlay_edge(utils::Hmesh::Edge *e, Eigen::RowVector2d shift,
                       Tangram &tangram) {
  e = e->next();
  while (!tangram.is_underlay_edge(e->vi, e->next()->vi)) {
    shift += tangram.tangram_mesh.V.row(e->next()->vi) -
             tangram.tangram_mesh.V.row(e->twin()->vi);
    e = e->twin()->next();
  }
  return {e, shift};
}

static double mvc_coord(const Eigen::RowVector2d &p,
                        const Eigen::RowVector2d &vi,
                        const Eigen::RowVector2d &vprev,
                        const Eigen::RowVector2d &vnext) {
  double ang1 = std::acos((vprev - p).normalized().dot((vi - p).normalized()));
  double ang2 = std::acos((vnext - p).normalized().dot((vi - p).normalized()));
  return (std::tan(ang1 / 2) + std::tan(ang2 / 2)) / (vi - p).norm();
}

void Tangram::compute_pleats_mvc() {
  _index_to_pleat_mvc_index = std::vector<int>(_vertices_positions.rows(), -1);
  for (auto pleat : _pleat_indices) {
    // Get an underlay edge bounding the pleat face.
    auto [underlay_edge, shift] =
        get_bounding_underlay_edge(tangram_mesh.verts[pleat].edge()->fi, *this);
    auto e = underlay_edge;

    // Go over the underlay face and compute the MVC.
    auto periodic_basis = compute_periodic_shift(_vertices_positions);
    auto projection = periodic_basis.inverse();
    PleatMVC mvc{.pleat_index = pleat};
    double total_weight = 0;
    do {
      auto [next_e, next_shift] = get_next_underlay_edge(e, shift, *this);
      mvc.coords.push_back(
          {.index = e->next()->vi,
           .weight =
               mvc_coord(tangram_mesh.V.row(pleat),
                         tangram_mesh.V.row(e->next()->vi) + shift,
                         tangram_mesh.V.row(e->vi) + shift,
                         tangram_mesh.V.row(next_e->next()->vi) + next_shift),
           .shift = shift * projection});
      total_weight += mvc.coords.back().weight;
      e = next_e;
      shift = next_shift;
    } while (e != underlay_edge);
    for (auto &coord : mvc.coords) {
      coord.weight /= total_weight;
    }
    _pleat_mvc.push_back(mvc);
    _index_to_pleat_mvc_index[pleat] = _pleat_mvc.size() - 1;
  }
  std::cout << "Computed " << _pleat_mvc.size() << " pleats\n";
}

utils::Hmesh Tangram::get_tangram_face(int start_triangle) {
  // Search for an underlay edge to start from.
  auto [underlay_edge, shift] =
      get_bounding_underlay_edge(start_triangle, *this);

  assert(underlay_edge && "No underlay edge found, should not be possible...");

  auto e = underlay_edge;
  std::vector<Eigen::Vector2d> new_verts;
  auto periodic_basis = compute_periodic_shift(_vertices_positions);
  do {
    auto [next_e, next_shift] = get_next_underlay_edge(e, shift, *this);
    shift = shift * periodic_basis.inverse();
    new_verts.push_back(tangram_mesh.V.row(e->vi) + shift * periodic_basis);
    new_verts.push_back(tangram_mesh.V.row(e->next()->vi) +
                        shift * periodic_basis);
    e = next_e;
    shift = next_shift;
  } while (e != underlay_edge);
  Eigen::MatrixXd new_V = convert::to_eig_mat(new_verts);
  Eigen::Vector2d mean = new_V.colwise().mean();
  new_verts.push_back(mean);
  std::vector<std::vector<int>> new_faces;
  for (int i = 0; i < new_verts.size() / 2; i++) {
    new_faces.push_back({i * 2, i * 2 + 1, (int)new_verts.size() - 1});
  }
  return utils::Hmesh(convert::to_eig_mat(new_verts), new_faces);
}

void Tangram::classify_triangles() {
  _face_classification =
      std::vector<FaceClassification>(_triangulated_faces.size(), NOT_SET);

  // Find the grid faces composing the underlay face.
  auto collect_faces = [&](int start_face) {
    std::set<int> faces;
    std::vector<int> stack = {start_face};
    while (!stack.empty()) {
      int face = stack.back();
      stack.pop_back();
      if (!faces.insert(face).second)
        continue;
      auto e = tangram_mesh.face(face).edge();
      for (int i = 0; i < 3; i++, e = e->next()) {
        if (e->twin() && !is_underlay_edge(e->vi, e->next()->vi)) {
          stack.push_back(e->twin()->fi);
        }
      }
    }
    return faces;
  };

  auto classify_faces = [&](const std::set<int> &faces) {
    std::set<int> seen_sl, seen_vertices;
    for (int face : faces) {
      for (int i = 0; i < 3; i++) {
        int v = _triangulated_faces[face][i];
        if (!seen_vertices.insert(v).second)
          continue;
        int sl = index_to_sl(v);
        if (sl != -1 && !seen_sl.insert(sl).second) {
          // Two of the underlay vertices are from the same stitching line.
          return PLEAT;
        }
      }
    }
    return UNDERLAY;
  };

  for (int i = 0; i < _triangulated_faces.size(); i++) {
    if (_face_classification[i] != NOT_SET)
      continue;
    auto faces = collect_faces(i);
    FaceClassification classification = classify_faces(faces);
    for (int face : faces) {
      _face_classification[face] = classification;
    }
  }
}

Eigen::MatrixXd Tangram::compute_vertices_positions(double Gary) {
  if (abs(Gary - 1.0) < 1e-6) {
    return _pattern.verts();
  }
  Optiz::Problem prob(_pattern.verts(),
                      {.cache_pattern = true, .remove_unreferenced = true});
  // Underlay edges should remain rigid.
  prob.add_element_energy(_underlay_edges.size(), [&](int i, auto &x) {
    int v0 = _underlay_edges[i][0];
    int v1 = _underlay_edges[i][1];
    double original_len =
        (_pattern.verts().row(v0) - _pattern.verts().row(v1)).squaredNorm();
    return Optiz::sqr((x.row(v1) - x.row(v0)).squaredNorm() - original_len);
  });

  // Stitching edges should shrink.
  auto stitching_edges = _pattern.get_stitching_pairs();
  prob.add_element_energy(stitching_edges.size(), [&](int i, auto &x) {
    int v0 = stitching_edges[i][0];
    int v1 = stitching_edges[i][1];
    double original_len =
        (_pattern.verts().row(v0) - _pattern.verts().row(v1)).squaredNorm();
    return Optiz::sqr((x.row(v1) - x.row(v0)).squaredNorm() -
                      original_len * Gary);
  });

  // TODO: Optimize for (nx() - 1) * (ny() - 1) vertices plus 2 translation
  // Align the y direction.
  prob.add_element_energy(_pattern.nx(), [&](int i, auto &x) {
    using T = FACTORY_TYPE(x);
    if (_coord_to_sl(i, 0) == -1) {
      return T(0.0);
    }
    int v0 = _pattern.coords_to_index(i, 0);
    int v1 = _pattern.coords_to_index(i, _pattern.ny() - 1);
    return Optiz::sqr(x(v0, 0) - x(v1, 0));
  });
  prob.add_element_energy(_pattern.ny(), [&](int i, auto &x) {
    using T = FACTORY_TYPE(x);
    if (_coord_to_sl(0, i) == -1) {
      return T(0.0);
    }
    int v0 = _pattern.coords_to_index(0, i);
    int v1 = _pattern.coords_to_index(_pattern.nx() - 1,i);
    return Optiz::sqr(x(v0, 1) - x(v1, 1));
  });

  // prob.options().set_report_level(Optiz::Problem::Options::NONE);

  Eigen::MatrixXd sol = prob.optimize().x();
  auto periodic_basis = compute_periodic_shift(sol);
  for (auto &pleat : _pleat_mvc) {
    Eigen::RowVector2d pos = Eigen::RowVector2d::Zero();
    for (const auto &coord : pleat.coords) {
      pos +=
          coord.weight * (sol.row(coord.index) + coord.shift * periodic_basis);
    }
    sol.row(pleat.pleat_index) = pos;
  }
  return sol;
}

void Tangram::set_stitching_scale(double stitching_scale) {
  _stitching_scale = stitching_scale;
  _vertices_positions = compute_vertices_positions(_stitching_scale);
}

Eigen::Matrix2d Tangram::compute_periodic_basis() const {
  return compute_periodic_shift(_vertices_positions);
}

Eigen::Matrix2d Tangram::compute_periodic_shift(const Eigen::MatrixXd& verts_positions) const {
  Eigen::Matrix2d res;
  for (int i = 0; i < _pattern.nx(); i++) {
    if (_coord_to_sl(i, 0) == -1) {
      continue;
    }
    res.row(1) = verts_positions.row(
                     _pattern.coords_to_index(i, _pattern.ny() - 1)) -
                 verts_positions.row(_pattern.coords_to_index(i, 0));
    break;
  }
  for (int i = 0; i < _pattern.ny(); i++) {
    if (_coord_to_sl(0, i) == -1) {
      continue;
    }
    res.row(0) = verts_positions.row(
                     _pattern.coords_to_index(_pattern.nx() - 1, i)) -
                 verts_positions.row(_pattern.coords_to_index(0, i));
    break;
  }
  return res;
}

} // namespace smocking