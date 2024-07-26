
#include "pullback.h"
#include "smocking_pattern/Tangram.h"
#include <igl/per_vertex_normals.h>
#include <polyscope/surface_mesh.h>
#include <parameterization/remeshing.h>
#include <state/state.h>
#include <utils/point_in_param.h>

namespace param {

Eigen::MatrixXd PullbackResult::get_open_configuration_no_symmetry(double opening) const {
  Eigen::MatrixXd verts(vertices.size(), 2);
  Eigen::MatrixXd tangram_verts =
      state::tangram.compute_vertices_positions(opening);
  Eigen::Matrix2d periodic_basis =
      state::tangram.compute_periodic_shift(tangram_verts);

  for (int i = 0; i < vertices.size(); i++) {
    if (!vertices[i].is_inside)
      continue;
    verts.row(i) = tangram_verts.row(vertices[i].unit_pattern_index) +
                   vertices[i].repeat * periodic_basis;
  }
  return verts;
}

Eigen::MatrixXd PullbackResult::get_closed_configuration_no_symmetry() const {
  auto& mesh = state::mesh;
  Eigen::MatrixXd verts(vertices.size(), 3);
  for (int i = 0; i < vertices.size(); i++) {
    if (!vertices[i].is_inside)
      continue;
    verts.row(i) = vertices[i].bary_coords.transpose() *
                   mesh.V(mesh.F[(vertices[i].uv_face)], Eigen::all);
  }
  // Get average underlay edge size for the scaling.
  double avg_underlay_edge_size = 0.0;
  for (auto &edge : underlay_edges) {
    avg_underlay_edge_size += (verts.row(edge[0]) - verts.row(edge[1])).norm();
  }
  avg_underlay_edge_size /= underlay_edges.size();
  for (int i = 0; i < vertices.size(); i++) {
    if (vertices[i].is_underlay || vertices[i].uv_face == -1)
      continue;
    Eigen::RowVector3d nrm =
        (vertices[i].bary_coords.transpose() *
         (mesh.Vn(mesh.F[(vertices[i].uv_face)], Eigen::all)))
            .normalized();
    verts.row(i) += nrm * avg_underlay_edge_size * 0.5;
  }
  return verts;
}

int local_index_to_global(int local_index, int i, int j, int global_nx) {
  auto &pattern = state::tangram.get_pattern();
  Eigen::Vector2i local_coords = pattern.index_to_coords(local_index);
  Eigen::Vector2i global_coords =
      local_coords +
      Eigen::Vector2i((pattern.nx() - 1) * i, (pattern.ny() - 1) * j);
  return global_coords.x() + global_coords.y() * global_nx;
}

static double mvc_coord(const Eigen::Vector2d &p, const Eigen::Vector2d &vi,
                        const Eigen::Vector2d &vprev,
                        const Eigen::Vector2d &vnext) {
  double ang1 = std::acos((vprev - p).normalized().dot((vi - p).normalized()));
  double ang2 = std::acos((vnext - p).normalized().dot((vi - p).normalized()));
  return (std::tan(ang1 / 2) + std::tan(ang2 / 2)) / (vi - p).norm();
}

void fill_mvc(const Eigen::Vector2d &vp, const std::vector<int> &inds,
              const PullbackResult &result, smocking::Tangram::PleatMVC &mvc) {
  for (int i = 0; i < inds.size(); i++) {
    Eigen::Vector2d vi = result.vertices[inds[i]].open_pos;
    Eigen::Vector2d vprev =
        result.vertices[inds[(i - 1 + inds.size()) % inds.size()]].open_pos;
    Eigen::Vector2d vnext =
        result.vertices[inds[(i + 1) % inds.size()]].open_pos;
    mvc.coords.push_back(
        {.index = inds[i], .weight = mvc_coord(vp, vi, vprev, vnext)});
  }
  // Normalize the weights.
  double total_weight = 0.0;
  for (auto &coord : mvc.coords) {
    total_weight += coord.weight;
  }
  for (auto &coord : mvc.coords) {
    coord.weight /= total_weight;
  }
}

void update_repeat_pleat_mvc(const Eigen::Vector2i &min_uv,
                             const Eigen::Vector2i &num_repeats, int i, int j,
                             PullbackResult &result) {
  auto &pattern = state::tangram.get_pattern();
  int nx = state::tangram.get_pattern().nx(),
      ny = state::tangram.get_pattern().ny();
  // Update vertices positions.
  int ii = min_uv(0) + i;
  int jj = min_uv(1) + j;
  int global_nx = (num_repeats(0) * (nx - 1) + 1);
  int global_ny = (num_repeats(1) * (ny - 1) + 1);
  for (auto &mvc : state::tangram.get_pleat_mvc()) {
    // Calculate the global index of the pleat.
    int pleat_index = local_index_to_global(mvc.pleat_index, i, j, global_nx);
    // Check if all the surrounding vertices are inside the UV.
    bool all_inside = true;
    std::vector<int> inside_polygon;
    smocking::Tangram::PleatMVC new_mvc{.pleat_index = pleat_index};
    for (auto &coord : mvc.coords) {
      int i_shift = coord.shift.x(), j_shift = coord.shift.y();
      if (i + i_shift < 0 || i + i_shift >= num_repeats(0) || j + j_shift < 0 ||
          j + j_shift >= num_repeats(1)) {
        all_inside = false;
        // break;
        continue;
      }
      int global_index = local_index_to_global(coord.index, i + i_shift,
                                               j + j_shift, global_nx);
      if (!result.vertices[global_index].is_inside) {
        all_inside = false;
        continue;
      }
      inside_polygon.push_back(global_index);
      // new_mvc.coords.push_back({global_index, coord.weight});
    }
    // TODO: Remove inside_polygon.size() > 2 check, only keep non-boundary pleats.
    if (all_inside /*|| inside_polygon.size() > 2*/) {
      fill_mvc(result.vertices[pleat_index].open_pos, inside_polygon, result,
               new_mvc);
      result.pleat_mvcs.push_back(new_mvc);
      result.vertices[pleat_index].is_boundary_pleat = !all_inside;
    } else {
      result.vertices[pleat_index].is_inside = false;
    }
  }
}

void update_faces(int i, int j, int global_nx,
                  const std::vector<int> orig_to_inside,
                  PullbackResult &result) {
  // Update faces inside the UV.
  for (int f = 0; f < state::tangram.get_triangles().size(); f++) {
    auto &face = state::tangram.get_triangles()[f];
    int v0 = local_index_to_global(face[0], i, j, global_nx),
        v1 = local_index_to_global(face[1], i, j, global_nx),
        v2 = local_index_to_global(face[2], i, j, global_nx);
    if (result.vertices[v0].is_inside && result.vertices[v1].is_inside &&
        result.vertices[v2].is_inside) {
      result.faces.push_back({v0, v1, v2});
      result.reduced_faces.push_back(
          {orig_to_inside[v0], orig_to_inside[v1], orig_to_inside[v2]});
      result.face_classification.push_back(
          state::tangram.get_face_classification()[f]);
    }
  }
}

void process_repeat(const Eigen::Vector2i &min_uv,
                    const Eigen::Vector2i &num_repeats, int i, int j,
                    const Eigen::Matrix2d &periodic_basis,
                    const Eigen::MatrixXd &tangram_verts,
                    const Eigen::MatrixXd &UV,
                    const std::vector<Eigen::Vector2i> &stitching_pairs,
                    PullbackResult &result) {
  auto &pattern = state::tangram.get_pattern();
  int nx = state::tangram.get_pattern().nx(),
      ny = state::tangram.get_pattern().ny();
  // Update vertices positions.
  int ii = min_uv(0) + i;
  int jj = min_uv(1) + j;
  int global_nx = (num_repeats(0) * (nx - 1) + 1);
  Eigen::RowVector2d shift = Eigen::RowVector2d(ii, jj) * periodic_basis;
  // Whether the vertices are inside the UV.
  Eigen::VectorX<bool> inside(nx * ny);
  for (int k = 0; k < tangram_verts.rows(); k++) {
    Eigen::Vector2d p = tangram_verts.row(k) + shift;
    int global_index = local_index_to_global(k, i, j, global_nx);
    // TODO: Can improve performance.
    auto res = utils::point_in_param(UV, state::mesh.F, p);
    inside(k) = res.first != -1;
    result.vertices[global_index] = {
        k,         {ii, jj},   p,
        res.first, res.second, state::tangram.is_underlay_vertex(k),
        inside(k)};
  }
  // Update underlay edges inside the UV.
  for (auto &edge : state::tangram.get_underlay_edges()) {
    int index1 = edge(0), index2 = edge(1);
    if (inside(index1) && inside(index2)) {
      result.underlay_edges.push_back(
          {local_index_to_global(index1, i, j, global_nx),
           local_index_to_global(index2, i, j, global_nx)});
    }
  }
  // Update stitching edges.
  for (auto &edge : stitching_pairs) {
    int index1 = edge(0), index2 = edge(1);
    if (inside(index1) && inside(index2)) {
      result.stitching_edges.push_back(
          {local_index_to_global(index1, i, j, global_nx),
           local_index_to_global(index2, i, j, global_nx)});
    }
  }
}

void pullback_no_rosy(double rotate_uv_angle, double scale,
                      const Eigen::Vector2d &shift) {
  assert(state::field_data.field.N == 2 && "Only implemented for N=2");
  // Close the tangram.
  state::tangram.set_stitching_scale(0);
  // Transform the UV coordinates.
  Eigen::MatrixXd UV = state::field_data.NFunction.leftCols(2);
  UV.rowwise() -= UV.colwise().mean();
  UV = UV * Eigen::Rotation2Dd(-rotate_uv_angle).toRotationMatrix();
  UV *= scale;
  UV.rowwise() += shift.transpose();
  // auto res = polyscope::registerSurfaceMesh2D("uv", UV, state::mesh.F);
  // auto para = res->addVertexParameterizationQuantity("uvuv", UV);
  // para->setEnabled(true);
  // para->setCheckerSize(1.0);

  // Calculate the bounds for the repetition of the unit tile.
  Eigen::Matrix2d periodic_basis = state::tangram.compute_periodic_basis();
  Eigen::MatrixXd transformed_uv = UV * periodic_basis.inverse();
  Eigen::Vector2i min_uv =
      transformed_uv.colwise().minCoeff().array().floor().cast<int>() - 1;
  Eigen::Vector2i max_uv =
      transformed_uv.colwise().maxCoeff().array().ceil().cast<int>() + 1;
  Eigen::Vector2i num_repeats = max_uv - min_uv;

  auto tangram_verts = state::tangram._vertices_positions;
  auto &pattern = state::tangram.get_pattern();
  auto stitching_pairs = pattern.get_stitching_pairs();
  int global_nx = (num_repeats(0) * (pattern.nx() - 1) + 1);
  int global_ny = (num_repeats(1) * (pattern.ny() - 1) + 1);
  Eigen::MatrixXd verts(global_nx * global_ny, 2);
  PullbackResult &result = state::pullback_result;
  result = std::move(PullbackResult());
  result.get_open_configuration = [&](double opening) {
    return result.get_open_configuration_no_symmetry(opening);
  };
  result.get_closed_configuration = [&]() {
    return result.get_closed_configuration_no_symmetry();
  };
  
  result.vertices.resize(global_nx * global_ny);
  for (int i = 0; i < num_repeats(0); i++) {
    for (int j = 0; j < num_repeats(1); j++) {
      process_repeat(min_uv, num_repeats, i, j, periodic_basis, tangram_verts,
                     UV, stitching_pairs, result);
    }
  }
  for (int i = 0; i < num_repeats(0); i++) {
    for (int j = 0; j < num_repeats(1); j++) {
      update_repeat_pleat_mvc(min_uv, num_repeats, i, j, result);
    }
  }
  // Indices of the vertices inside the UV.
  result.orig_to_inside = std::vector<int>(result.vertices.size());
  for (int v = 0; v < result.vertices.size(); v++) {
    if (result.vertices[v].is_inside) {
      result.inside_verts_indices.push_back(v);
      if (result.vertices[v].is_underlay) {
        result.inside_underlay_indices.push_back(
            result.inside_verts_indices.size() - 1);
      }
      result.orig_to_inside[v] = result.inside_verts_indices.size() - 1;
    }
  }
  for (int i = 0; i < num_repeats(0); i++) {
    for (int j = 0; j < num_repeats(1); j++) {
      update_faces(i, j, global_nx, result.orig_to_inside, result);
    }
  }
}

} // namespace param