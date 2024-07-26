#include "optimize.h"
#include "Optiz/NewtonSolver/Problem.h"
#include "utils/Hmesh.h"
#include <Optiz/Problem.h>
#include <polyscope/curve_network.h>
#include <state/state.h>
#include <vector>

namespace layout {

bool stop_optimization = true;
bool layout_updated;

struct SliceIndices {
  std::vector<int> new_to_old;
  std::vector<int> old_to_new;
};
SliceIndices get_underlay_indices() {
  SliceIndices underlay_indices;
  underlay_indices.old_to_new.resize(state::pullback_result.vertices.size(),
                                     -1);
  for (int i = 0; i < state::pullback_result.vertices.size(); i++) {
    if (state::pullback_result.vertices[i].is_underlay &&
        state::pullback_result.vertices[i].uv_face != -1) {
      underlay_indices.old_to_new[i] = underlay_indices.new_to_old.size();
      underlay_indices.new_to_old.push_back(i);
    }
  }
  return underlay_indices;
}

template <typename T>
T cross2d(const Eigen::VectorX<T> &v1, const Eigen::VectorX<T> &v2) {
  return v1.x() * v2.y() - v1.y() * v2.x();
}
template <typename T>
T angle_between(const Eigen::VectorX<T> &v1, const Eigen::VectorX<T> &v2) {
  return atan2(cross2d(v1, v2), v1.dot(v2));
}
struct PleatConstraint {
  // Angle between (v1 - v0), (v2 - v1).
  int v0, v1, v2;
  double angle;
};

std::vector<PleatConstraint>
get_pleat_constraints(const Eigen::MatrixXd &init) {

  std::vector<PleatConstraint> constraints;
  for (auto &pleat_mvc : state::pullback_result.pleat_mvcs) {
    if (pleat_mvc.coords.size() < 3)
      continue;
    int is_open = state::pullback_result.vertices[pleat_mvc.pleat_index]
                      .is_boundary_pleat;
    int limit = is_open ? pleat_mvc.coords.size() - 2 : pleat_mvc.coords.size();
    for (int i = 0; i < limit; i++) {
      int v0 = pleat_mvc.coords[i].index;
      int v1 = pleat_mvc.coords[(i + 1) % pleat_mvc.coords.size()].index;
      int v2 = pleat_mvc.coords[(i + 2) % pleat_mvc.coords.size()].index;
      Eigen::VectorXd e10 = (init.row(v1) - init.row(v0)).normalized(),
                      e21 = (init.row(v2) - init.row(v1)).normalized();
      constraints.push_back(PleatConstraint{
          .v0 = v0, .v1 = v1, .v2 = v2, .angle = angle_between(e10, e21)});
    }
  }
  return constraints;
}

struct SeamlessConstraints {
  // (x(v0) - x(v1)).norm() == (x(tv0) - x(tv1)).norm().
  struct EqualLengthConstraint {
    int v0, v1;
    int tv0, tv1;
  };
  // Angle between (v1 - v0), (v2 - v1) == angle between (tv1 - tv0), (tv2 -
  // tv1).
  struct AngleConstraint {
    int v0, v1, v2;
    int tv0, tv1, tv2;
  };
  std::vector<EqualLengthConstraint> equal_length_constraints;
  std::vector<AngleConstraint> angle_constraints;
};
SeamlessConstraints get_seamless_constraints(const Eigen::MatrixXd &verts) {
  SeamlessConstraints result;
  if (state::pullback_result.duplicate_verts.empty())
    return result;
  utils::Hmesh mesh(verts, state::pullback_result.faces);
  auto &dup_verts = state::pullback_result.duplicate_verts;
  // Returns the next boundary vertex.
  auto next_bv = [&](int v) { return mesh.verts[v].edge()->next()->vi; };
  auto is_straight = [&](int v0, int v1, int v2) {
    Eigen::VectorXd e1 = verts.row(v1) - verts.row(v0),
                    e2 = verts.row(v2) - verts.row(v1);
    return abs(angle_between(e1, e2)) < 1e-3;
  };
  auto dup = [&](int v) { return *dup_verts[v].begin(); };

  // Iterate over underlay duplicate vertices.
  std::set<int> covered_len_verts;
  for (auto v : dup_verts) {
    int v0 = v.first, tv0 = dup(v0);
    if (!state::pullback_result.vertices[v0].is_underlay)
      continue;
    int v1 = next_bv(v0), v2 = next_bv(v1);
    if (!dup_verts.count(v2) || !is_straight(v0, v1, v2))
      continue;
    if (covered_len_verts.insert(v0).second) {
      result.equal_length_constraints.emplace_back(v0, v2, dup(v0), dup(v2));
      covered_len_verts.insert(dup(v2));
    } else {
      continue;
    }

    // Angle constraints.
    int v3 = next_bv(v2), v4 = next_bv(v3);
    if (!dup_verts.count(v4) || !is_straight(v2, v3, v4))
      continue;
    result.angle_constraints.emplace_back(v0, v2, v4, dup(v0), dup(v2),
                                          dup(v4));
  }

  return result;
}

static void fill_pleats(Eigen::MatrixXd &layout) {
  for (auto &pleat_mvc : state::pullback_result.pleat_mvcs) {
    Eigen::RowVector2d pleat_pos = Eigen::RowVector2d::Zero();
    for (auto &coord : pleat_mvc.coords) {
      pleat_pos += layout.row(coord.index) * coord.weight;
    }
    layout.row(pleat_mvc.pleat_index) = pleat_pos;
  }
}

void optimize_for_surface() {
  Eigen::MatrixXd layout = state::pullback_result.get_open_configuration(1.0);
  Eigen::MatrixXd surface_verts =
      state::pullback_result.get_closed_configuration();
  // Scale the layout to match the surface.
  double avg_scale_diff = 0.0;
  for (auto &edge : state::pullback_result.underlay_edges) {
    avg_scale_diff +=
        (layout.row(edge[0]) - layout.row(edge[1])).norm() /
        (surface_verts.row(edge[0]) - surface_verts.row(edge[1])).norm();
  }
  avg_scale_diff /= state::pullback_result.underlay_edges.size();
  layout /= avg_scale_diff;

  // Initialize pleat_weight to 100.
  double pleat_weight = 100;

  // Only optimize the underlay vertices.
  auto udnerlay_indices = get_underlay_indices();
  auto new_ind = [&](int i) { return udnerlay_indices.old_to_new[i]; };
  Optiz::Problem prob(layout(udnerlay_indices.new_to_old, Eigen::all));

  // Underlay energy.
  auto edge_len_energy = [&](int i, auto &x) {
    auto edge = state::pullback_result.underlay_edges[i];
    int v0 = edge(0), v1 = edge(1);
    double surface_edge_len =
        (surface_verts.row(v0) - surface_verts.row(v1)).norm();
    auto current_length = (x.row(new_ind(v0)) - x.row(new_ind(v1))).norm();
    return Optiz::sqr(current_length / surface_edge_len - 1);
  };
  auto total_edge_energy = [&] {
    double res = 0.0;
    auto x = prob.x();
    for (int i = 0; i < state::pullback_result.underlay_edges.size(); i++)
      res += edge_len_energy(i, x);
    return res;
  };

  // Add the energy term for the underlay edges.
  prob.add_element_energy<4>(state::pullback_result.underlay_edges.size(),
                             edge_len_energy);

  auto pleat_constraints = get_pleat_constraints(layout);
  prob.add_element_energy<6>(pleat_constraints.size(), [&](int i, auto &x) {
    using T = FACTORY_TYPE(x);
    auto &constraint = pleat_constraints[i];
    auto x0 = x.row(new_ind(constraint.v0)), x1 = x.row(new_ind(constraint.v1)),
         x2 = x.row(new_ind(constraint.v2));
    Eigen::VectorX<T> e10 = x1 - x0, e21 = x2 - x1;
    return pleat_weight *
           Optiz::sqr((angle_between(e10, e21) - constraint.angle) /
                      (2 * M_PI));
  });

  auto seam_constraints = get_seamless_constraints(layout);
  if (!seam_constraints.equal_length_constraints.empty()) {
    prob.add_element_energy<8>(
        seam_constraints.equal_length_constraints.size(), [&](int i, auto &x) {
          auto &constraint = seam_constraints.equal_length_constraints[i];
          auto x0 = x.row(new_ind(constraint.v0)),
               x1 = x.row(new_ind(constraint.v1));
          auto y0 = x.row(new_ind(constraint.tv0)),
               y1 = x.row(new_ind(constraint.tv1));
          auto x01_norm = (x1 - x0).norm();
          auto y01_norm = (y1 - y0).norm();
          return 0.1 * Optiz::sqr(x01_norm / y01_norm - 1);
        });
  }
  if (!seam_constraints.angle_constraints.empty()) {
    prob.add_element_energy<12>(
        seam_constraints.angle_constraints.size(), [&](int i, auto &x) {
          using T = FACTORY_TYPE(x);
          auto &constraint = seam_constraints.angle_constraints[i];
          auto x0 = x.row(new_ind(constraint.v0)),
               x1 = x.row(new_ind(constraint.v1)),
               x2 = x.row(new_ind(constraint.v2));
          auto y0 = x.row(new_ind(constraint.tv0)),
               y1 = x.row(new_ind(constraint.tv1)),
               y2 = x.row(new_ind(constraint.tv2));
          auto angle1 = angle_between<T>((x1 - x0), (x2 - x1)),
               angle2 = angle_between<T>((y1 - y0), (y2 - y1));
          return 0.1 * Optiz::sqr((angle1 - angle2) / (2 * M_PI));
        });
  }
  // Optimize and update the underlay vertices in the layout.
  prob.options().remove_unreferenced = true;
  prob.options().cache_pattern = true;
  prob.options().set_report_level(Optiz::Problem::Options::NONE);

  // std::vector<Eigen::Vector2i> edges;
  // std::vector<int> quantity;
  // for (auto &bla : seam_constraints.equal_length_constraints) {
  //   edges.push_back({bla.v0, bla.v1});
  //   edges.push_back({bla.tv0, bla.tv1});
  //   int cur = quantity.size();
  //   quantity.push_back(cur);
  //   quantity.push_back(cur);
  // }
  // auto e = polyscope::registerCurveNetwork2D("bla", layout, edges);
  // e->addEdgeScalarQuantity("quant", quantity);
  // return;

  stop_optimization = false;
  for (int i = 0; i < 70 && !stop_optimization; i++) {
    prob.optimize();
    layout(udnerlay_indices.new_to_old, Eigen::all) = prob.x();
    // Update pleats with mean value coordinates.
    fill_pleats(layout);
    state::optimized_layout = layout;
    double current_underlay_energy = total_edge_energy();
    std::cout << "Iteration " << i << ": " << current_underlay_energy
              << std::endl;
    pleat_weight *= 0.6;
    layout_updated = true;
    if (current_underlay_energy < 1e-4) {
      break;
    }
  }
  stop_optimization = true;
}
} // namespace layout