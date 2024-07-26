#include "smocking.h"
#include "manager.h"
#include "polyscope/group.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/types.h"
#include "utils/conversions.h"
#include <polyscope/curve_network.h>
#include <state/state.h>
#include <utils/replicate.h>

namespace view {

const Eigen::Vector3d underlay_region_color = {0.2823, 0.6274, 0.8353};
// const Eigen::Vector3d underlay_region_color = {195/255.0, 226/255.0,
// 245/255.0};
const Eigen::Vector3d underlay_edges_color = {133 / 255.0, 198 / 255.0,
                                              237 / 255.0};
const Eigen::Vector3d pleat_region_color = {0.8980, 0.1922, 0.4039};
// const Eigen::Vector3d pleat_region_color = {254/255.0, 212/255.0, 228/255.0};
const Eigen::Vector3d pleat_edges_color = {99 / 255.0, 100 / 255.0,
                                           102 / 255.0};

void visualize_smocking_pattern(const smocking::UnitSmockingPattern &pattern) {
  hide_all();
  auto group = polyscope::createGroup("smocking");
  auto network = polyscope::registerCurveNetwork2D(
      "smocking pattern", pattern.verts(), pattern.get_stitching_pairs());
  group->addChildStructure(*network);
}

std::tuple<Eigen::MatrixXd, std::vector<Eigen::Vector2i>>
remove_pleat_verts(const Eigen::MatrixXd &verts,
                   const std::vector<Eigen::Vector2i> &edges) {
  std::vector<int> keep_verts;
  std::vector<int> old_to_new(verts.rows(), -1);
  std::vector<Eigen::Vector2i> new_edges;
  for (auto &edge : edges) {
    if (old_to_new[edge[0]] == -1) {
      old_to_new[edge[0]] = keep_verts.size();
      keep_verts.push_back(edge[0]);
    }
    if (old_to_new[edge[1]] == -1) {
      old_to_new[edge[1]] = keep_verts.size();
      keep_verts.push_back(edge[1]);
    }
    new_edges.push_back({old_to_new[edge[0]], old_to_new[edge[1]]});
  }
  return {verts(keep_verts, Eigen::all).eval(), new_edges};
}

void visualize_tangram(const smocking::Tangram &tangram, int repeat_x,
                       int repeat_y) {
  hide_all();

  Eigen::MatrixXd verts = convert::to_3d(tangram._vertices_positions);
  for (auto pleat : tangram._pleat_indices) {
    verts(pleat, 2) = -1.0 + tangram._stitching_scale;
  }

  auto periodic_basis = tangram.compute_periodic_basis();
  // Get all the vertices positions.
  Eigen::MatrixXd all_verts = verts.replicate(repeat_x * repeat_y, 1);
  for (int x = 0; x < repeat_x; x++) {
    for (int y = 0; y < repeat_y; y++) {
      int ind = x * repeat_y + y;
      all_verts.block(ind * verts.rows(), 0, verts.rows(), 2).rowwise() +=
          Eigen::RowVector2d(x, y) * periodic_basis;
    }
  }
  all_verts.rowwise() -= all_verts.colwise().mean();

  polyscope::Group *group = polyscope::state::groups.count("smocking") == 0
                                ? polyscope::createGroup("smocking")
                                : polyscope::getGroup("smocking");
  // Stitching lines.
  // Get all stitching pairs.
  auto all_stitching_pairs =
      replicate::vec_with_shift(tangram.get_pattern().get_stitching_pairs(),
                                repeat_x * repeat_y, (int)verts.rows());
  auto [keep_verts, keep_pairs] =
      remove_pleat_verts(all_verts, all_stitching_pairs);
  auto stitch_edges = polyscope::registerCurveNetwork("stitching lines",
                                                      keep_verts, keep_pairs);
  stitch_edges->setRadius(0.015, false);
  stitch_edges->setColor(convert::to_glm(pleat_edges_color));
  // stitch_edges->setMaterial("flat");
  group->addChildStructure(*stitch_edges);
  // Underlay edges.
  auto all_underlay_edges = replicate::vec_with_shift(
      tangram.get_underlay_edges(), repeat_x * repeat_y, (int)verts.rows());
  auto [keep_verts2, keep_underlay] =
      remove_pleat_verts(all_verts, all_underlay_edges);
  auto underlay_edges = polyscope::registerCurveNetwork(
      "underlay edges", keep_verts2, keep_underlay);
  underlay_edges->setRadius(0.015, false);
  underlay_edges->setColor(convert::to_glm(underlay_edges_color));
  // underlay_edges->setMaterial("flat");
  group->addChildStructure(*underlay_edges);

  // Tangram mesh.
  auto all_tangram_triangles = replicate::vec_with_shift(
      tangram.get_triangles(), repeat_x * repeat_y, (int)verts.rows());
  auto faces = polyscope::registerSurfaceMesh("tangram faces", all_verts,
                                              all_tangram_triangles);
  group->addChildStructure(*faces);
  // faces->setMaterial("flat");

  // Color its faces.
  Eigen::MatrixXd mesh_colors =
      Eigen::MatrixXd::Zero(all_tangram_triangles.size(), 3);
  for (int i = 0; i < mesh_colors.rows(); i++) {
    switch (
        tangram.get_face_classification()[i % tangram.get_triangles().size()]) {
    case smocking::Tangram::FaceClassification::UNDERLAY:
      mesh_colors.row(i) = underlay_region_color;
      break;
    case smocking::Tangram::FaceClassification::PLEAT:
      mesh_colors.row(i) = pleat_region_color;
      break;
    default:
      break;
    }
  }
  faces->addFaceColorQuantity("face classification", mesh_colors)
      ->setEnabled(true);

  group->setEnabled(true);
}

void show_flat_tangram(Eigen::MatrixXd verts, const std::string &name) {
  remove_all();
  auto inside_verts =
      verts(state::pullback_result.inside_verts_indices, Eigen::all);
  auto tangram_mesh = polyscope::registerSurfaceMesh(
      name, inside_verts, state::pullback_result.reduced_faces);
  tangram_mesh->setEnabled(true);
  Eigen::MatrixXd mesh_colors =
      Eigen::MatrixXd::Zero(state::pullback_result.faces.size(), 3);
  for (int i = 0; i < mesh_colors.rows(); i++) {
    switch (state::pullback_result.face_classification[i]) {
    case smocking::Tangram::FaceClassification::UNDERLAY:
      mesh_colors.row(i) = underlay_region_color;
      break;
    case smocking::Tangram::FaceClassification::PLEAT:
      mesh_colors.row(i) = pleat_region_color;
      break;
    default:
      break;
    }
  }
  tangram_mesh->addFaceColorQuantity("face classification", mesh_colors)
      ->setEnabled(true);
}

void show_pulled_back_tangram() {
  remove_all();
  Eigen::MatrixXd verts = state::pullback_result.get_closed_configuration();
  show_flat_tangram(verts, "closed configuration");
}

void show_open_configuration() {
  remove_all();
  Eigen::MatrixXd verts =
      convert::to_3d(state::pullback_result.get_open_configuration(1.0));
  Eigen::MatrixXd surface_verts =
      state::pullback_result.get_closed_configuration();
  // Scale the layout to match the surface.
  double avg_scale_diff = 0.0;
  for (auto &edge : state::pullback_result.underlay_edges) {
    avg_scale_diff +=
        (verts.row(edge[0]) - verts.row(edge[1])).norm() /
        (surface_verts.row(edge[0]) - surface_verts.row(edge[1])).norm();
  }
  avg_scale_diff /= state::pullback_result.underlay_edges.size();
  verts /= avg_scale_diff;
  show_flat_tangram(verts, "open configuration");
}

void show_optimized_layout() {
  remove_all();
  Eigen::MatrixXd verts = convert::to_3d(state::optimized_layout);
  show_flat_tangram(verts, "optimized layout");
  // hide_all();
  // auto tangram_mesh = polyscope::registerSurfaceMesh2D(
  //     "optimized_layout",
  //     state::optimized_layout(state::pullback_result.inside_verts_indices,
  //                             Eigen::all),
  //     state::pullback_result.reduced_faces);
  // tangram_mesh->setEnabled(true);
  // Eigen::MatrixXd mesh_colors =
  //     Eigen::MatrixXd::Zero(state::pullback_result.faces.size(), 3);
  // for (int i = 0; i < mesh_colors.rows(); i++) {
  //   switch (state::pullback_result.face_classification[i]) {
  //   case smocking::Tangram::FaceClassification::UNDERLAY:
  //     mesh_colors.row(i) = underlay_region_color;
  //     break;
  //   case smocking::Tangram::FaceClassification::PLEAT:
  //     mesh_colors.row(i) = pleat_region_color;
  //     break;
  //   default:
  //     break;
  //   }
  // }
  // tangram_mesh->addFaceColorQuantity("face classification", mesh_colors)
  //     ->setEnabled(true);
}

void show_stitching_lines() {
  Eigen::MatrixXd verts = convert::to_3d(state::optimized_layout);
  std::set<int> keep_verts;
  for (auto &face : state::pullback_result.faces) {
    for (auto v : face) {
      keep_verts.insert(v);
    }
  }
  std::vector<Eigen::Vector2i> edges;
  for (auto &edge : state::pullback_result.stitching_edges) {
    if (keep_verts.count(edge[0]) && keep_verts.count(edge[1])) {
      edges.push_back(Eigen::Vector2i(edge[0], edge[1]));
    }
  }
  auto [v, e] = remove_pleat_verts(verts, edges);
  auto stitch_edges = polyscope::registerCurveNetwork("stitching_lines", v, e);
  stitch_edges->setRadius(0.003, false);
  stitch_edges->setColor(glm::vec3(0, 0, 0));
}

void show_arap_result() {
  remove_all();
  auto arap_mesh = polyscope::registerSurfaceMesh("arap", state::arap_verts,
                                                  state::pleats_arap->faces());
  arap_mesh->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
  arap_mesh->setSurfaceColor(glm::vec3(114.0 / 255, 135.0 / 255, 217.0 / 255));
}

} // namespace view