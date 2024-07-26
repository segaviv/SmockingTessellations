#include "auxetic_pattern.h"
#include "parameterization/pullback.h"
#include "utils/Hmesh.h"
#include "utils/construct.h"
#include <cmath>
#include <queue>
#include <utils/conversions.h>

namespace cut {

double opening_ang = M_PI * 2 / 3;

static const utils::Hmesh::Edge *
get_non_boundary_edge(const utils::Hmesh &mesh) {
  for (const auto &e : mesh.edges) {
    if (!e.is_boundary()) {
      return &e;
    }
  }
  return nullptr;
}

Eigen::MatrixXd
get_open_confiugration(const utils::Hmesh &orig_mesh, int nv,
                       const std::vector<std::vector<int>> &faces,
                       const Eigen::MatrixXi &is_cut_edge,
                       const std::vector<int> &new_to_old,
                       const std::vector<int> &pleats) {
  // Initialize the first edge to (0, 0) -> (1, 0).
  Eigen::MatrixXd result = Eigen::MatrixXd::Constant(nv, 2, 0);
  result.row(faces[0][1]) = Eigen::Vector2d(1, 0);
  auto ei = orig_mesh.verts_edges[orig_mesh.F[0][0]].at(orig_mesh.F[0][1]);
  std::set<const utils::Hmesh::Edge *> visited;
  std::queue<const utils::Hmesh::Edge *> edges_q;
  edges_q.push(&orig_mesh.edges[ei]);
  std::set<int> covered_verts;
  covered_verts.insert(faces[0][0]);
  covered_verts.insert(faces[0][1]);
  auto update_result = [&](int v, const auto &x) {
    if (covered_verts.insert(v).second) {
      result.row(v) = x;
    }
  };

  // Sine law.
  double pleat_edge_length =
      std::sin(M_PI * 2 / 3 - opening_ang / 2) / std::sin(M_PI / 3);
  auto fv = [&](auto e) -> int { return faces[e->fi][e->fvi]; };
  auto pleat = [&](int v) -> int { return pleats[new_to_old[v]]; };

  Eigen::Rotation2Dd rot120(M_PI * 2 / 3);
  Eigen::Matrix2d opening_ang_mat =
      Eigen::Rotation2Dd(opening_ang).toRotationMatrix();
  Eigen::Matrix2d half_opening_ang_mat =
      Eigen::Rotation2Dd(opening_ang / 2).toRotationMatrix();
  while (!edges_q.empty()) {
    const utils::Hmesh::Edge *edge = edges_q.front();
    edges_q.pop();
    if (visited.insert(edge).second == false)
      continue;

    int v0 = fv(edge), v1 = fv(edge->next()), v2 = fv(edge->prev());
    Eigen::Vector2d edge_vec = (result.row(v1) - result.row(v0)).normalized();
    // Assign next vertex and push edge.
    update_result(v2, result.row(v1) + (rot120 * edge_vec).transpose());
    edges_q.push(edge->next());
    // Handle twin edge.
    if (is_cut_edge(edge->fi, edge->fvi)) {
      update_result(pleat(v0),
                    result.row(v1) +
                        (half_opening_ang_mat * (-edge_vec) * pleat_edge_length)
                            .transpose());
    } else {
      update_result(pleat(v1),
                    result.row(v0) + (half_opening_ang_mat.transpose() *
                                      edge_vec * pleat_edge_length)
                                         .transpose());
    }
    if (!edge->twin())
      continue;

    if ((!is_cut_edge(edge->fi, edge->fvi))) {
      update_result(fv(edge->twin()),
                    result.row(v0) +
                        (opening_ang_mat.transpose() * edge_vec).transpose());
    } else {
      update_result(fv(edge->twin()->next()),
                    result.row(v1) +
                        (opening_ang_mat * (-edge_vec)).transpose());
    }
    edges_q.push(edge->twin());
  }
  return result;
}

static double mvc_coord(const Eigen::RowVector2d &p,
                        const Eigen::RowVector2d &vi,
                        const Eigen::RowVector2d &vprev,
                        const Eigen::RowVector2d &vnext) {
  double ang1 = std::acos((vprev - p).normalized().dot((vi - p).normalized()));
  double ang2 = std::acos((vnext - p).normalized().dot((vi - p).normalized()));
  return (std::tan(ang1 / 2) + std::tan(ang2 / 2)) / (vi - p).norm();
}

void fill_mvc_coords(param::PullbackResult &result, int num_underlay_vertices,
                     int num_underlay_faces) {
  utils::Hmesh mesh(result.get_open_configuration(1.0), result.faces);
  for (int i = num_underlay_vertices; i < mesh.verts.size(); i++) {
    auto e = mesh.verts[i].edge();
    // Get underlay edge.
    while (!e->twin() || e->twin()->fi >= num_underlay_faces) {
      e = e->next();
    }
    auto e2 = e;
    if (e->prev()->twin()) {
      e2 = e->prev()->twin()->prev();
      while (e2 != e && e2->prev()->twin()) {
        e2 = e2->prev()->twin()->prev();
      }
    }
    if (!e2->prev()->twin()) {
      result.pleat_mvcs.push_back({.pleat_index = i});
      result.pleat_mvcs.back().coords.push_back(
          {.index = e2->vi, .weight = 0.5});
      while (e2->next()->twin()) {
        e2 = e2->next()->twin()->next();
        result.pleat_mvcs.back().coords.push_back(
            {.index = e2->vi, .weight = 0});
      }
      result.pleat_mvcs.back().coords.push_back(
          {.index = e2->next()->vi, .weight = 0.5});
      result.vertices[i].is_boundary_pleat = true;
      if (result.pleat_mvcs.back().coords.size() == 5) {
        result.pleat_mvcs.back().coords[0].weight = 0.25;
        result.pleat_mvcs.back().coords[1].weight = 0.25;
        result.pleat_mvcs.back().coords[3].weight = 0.25;
        result.pleat_mvcs.back().coords[4].weight = 0.25;
      }
      continue;
    }
    // Starting from e2, calculate mvc coordinates.
    result.pleat_mvcs.push_back({.pleat_index = i});
    do {
      int vprev = e->prev()->twin()->prev()->vi, vi = e->vi,
          vnext = e->next()->vi;
      result.pleat_mvcs.back().coords.push_back(
          {.index = vi,
           .weight = mvc_coord(mesh.V.row(i), mesh.V.row(vi), mesh.V.row(vprev),
                               mesh.V.row(vnext))});
      e = e->next()->twin()->next();
    } while (e != e2);
    result.pleat_mvcs.back().normalize_coords();
  }
}

void fill_stitch_verts(param::PullbackResult &result,
                       const utils::Hmesh &new_mesh,
                       const std::vector<std::set<int>> &old_to_new,
                       const std::vector<int> &new_to_old,
                       const std::vector<int> &pleats) {
  auto pleat = [&](int i) { return pleats[new_to_old[i]]; };
  auto edge_vi = [&](auto e) { return result.faces[e->fi][e->fvi]; };
  auto edge_tvi = [&](auto e) { return result.faces[e->fi][e->next()->fvi]; };
  utils::Hmesh mesh(result.get_closed_configuration(), result.faces);
  for (auto &edge : new_mesh.edges) {
    if (!edge.twin() || (edge.vi == edge.twin()->next()->vi &&
                         edge.next()->vi == edge.twin()->vi))
      continue;
    int vi = edge_vi(&edge), tvi = edge_tvi(edge.twin());
    if (vi != tvi) {
      if (mesh.is_boundary_vertex(vi) && mesh.is_boundary_vertex(tvi)) {
        result.duplicate_verts[vi].insert(tvi);
        result.duplicate_verts[tvi].insert(vi);
      }
      int pvi = pleat(vi), ptvi = pleat(tvi);
      if (pvi != ptvi) {
        result.duplicate_verts[pvi].insert(ptvi);
        result.duplicate_verts[ptvi].insert(pvi);
      }
    }
    vi = edge_tvi(&edge);
    tvi = edge_vi(edge.twin());
    if (vi != tvi) {
      if (mesh.is_boundary_vertex(vi) && mesh.is_boundary_vertex(tvi)) {
        result.duplicate_verts[vi].insert(tvi);
        result.duplicate_verts[tvi].insert(vi);
      }
      int pvi = pleat(vi), ptvi = pleat(tvi);
      if (pvi != ptvi) {
        result.duplicate_verts[pvi].insert(ptvi);
        result.duplicate_verts[ptvi].insert(pvi);
      }
    }
  }
  for (auto &e : result.duplicate_verts) {
    if (e.second.size() > 1) {
      std::cout << "Dup verts for " << e.first << ": ";
      for (auto v : e.second)
        std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "is it a pleat? "
                << (result.vertices[e.first].is_underlay ? "no" : "yes")
                << std::endl;
      std::cout << "Dups of first dup ( " << *e.second.begin() << " ): ";
      for (auto v : result.duplicate_verts[*e.second.begin()]) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
    }
  }
}

param::PullbackResult fill_pullback_result(
    const utils::Hmesh &orig_mesh, const utils::Hmesh &new_mesh,
    const std::vector<Eigen::VectorXd> &verts,
    const std::vector<std::vector<int>> &faces, int num_underlay,
    const std::vector<std::set<int>> &old_to_new,
    const std::vector<int> &new_to_old, const std::vector<int> &pleats,
    const std::set<const utils::Hmesh::Edge *> &covered_edges) {
  param::PullbackResult result;
  Eigen::MatrixXd V = convert::to_eig_mat(verts);
  int nv = (int)verts.size();
  result.get_closed_configuration = [V] { return V; };
  // Whether edge F[i][j] -> F[i][j + 1] is cut edge (F[i][j] is split).
  Eigen::MatrixXi is_cut_edge = Eigen::MatrixXi::Zero(faces.size(), 3);
  for (auto &e : covered_edges)
    is_cut_edge(e->fi, e->fvi) = 1;

  result.get_open_configuration = [orig_mesh, nv, faces, is_cut_edge,
                                   new_to_old, pleats](double t) {
    return get_open_confiugration(orig_mesh, nv, faces, is_cut_edge, new_to_old,
                                  pleats);
  };
  for (int i = 0; i < verts.size(); i++) {
    result.vertices.push_back(
        {.is_underlay = i < num_underlay, .is_inside = true});
  }
  result.faces = faces;

  // ARAP and visualization stuff.
  result.inside_verts_indices = construct::linspace((int)verts.size());
  result.inside_underlay_indices = construct::linspace(num_underlay);
  result.reduced_faces = result.faces;
  result.face_classification.resize(result.faces.size());
  for (int i = 0; i < result.faces.size(); i++) {
    result.face_classification[i] = i < orig_mesh.faces.size()
                                        ? param::FaceClassification::UNDERLAY
                                        : param::FaceClassification::PLEAT;
  }
  // Underlay edges.
  for (int i = 0; i < orig_mesh.faces.size(); i++) {
    for (int j = 0; j < 3; j++) {
      result.underlay_edges.emplace_back(result.faces[i][j],
                                         result.faces[i][(j + 1) % 3]);
    }
  }
  // Stitching edges.
  for (auto &s : old_to_new) {
    for (auto &v0 : s) {
      for (auto &v1 : s) {
        if (v0 > v1) {
          result.stitching_edges.emplace_back(v0, v1);
        }
      }
    }
  }
  fill_mvc_coords(result, num_underlay, orig_mesh.faces.size());
  fill_stitch_verts(result, new_mesh, old_to_new, new_to_old, pleats);
  return result;
}

static utils::Hmesh add_cut_twin_edges(const utils::Hmesh &mesh) {
  utils::Hmesh result(mesh);
  std::set<utils::Hmesh::Edge *> boundary_edges;
  for (auto &e : result.edges) {
    if (e.is_boundary())
      boundary_edges.insert(&e);
  }
  auto find_twin = [&](auto edge) -> utils::Hmesh::Edge * {
    for (auto e : boundary_edges) {
      if ((result.V.row(edge->vi) - result.V.row(e->next()->vi)).norm() <
              1e-6 &&
          (result.V.row(edge->next()->vi) - result.V.row(e->vi)).norm() <
              1e-6) {
        return e;
      }
    }
    return nullptr;
  };
  for (auto e : boundary_edges) {
    auto twin = find_twin(e);
    if (!twin || e->twin())
      continue;
    twin->ti = e->index;
    e->ti = twin->index;
  }
  return result;
}

bool has_real_twin(const utils::Hmesh::Edge *e) {
  return e->twin() && e->vi == e->twin()->next()->vi &&
         e->next()->vi == e->twin()->vi;
}

param::PullbackResult cut_auxetic_pattern(const utils::Hmesh &orig_mesh) {
  auto mesh = add_cut_twin_edges(orig_mesh);
  // BFS stuff.
  std::vector<std::vector<int>> faces = mesh.F;
  std::vector<Eigen::VectorXd> verts = convert::to_vec_mat(mesh.V);
  // Edges that have been visited.
  std::set<const utils::Hmesh::Edge *> covered_edges;
  std::queue<const utils::Hmesh::Edge *> edges_q;
  edges_q.push(get_non_boundary_edge(mesh));

  // Old to new and new to old mappings.
  std::vector<std::set<int>> old_to_new_vertices(mesh.nv());
  for (int i = 0; i < old_to_new_vertices.size(); i++) {
    old_to_new_vertices[i] = {i};
  }
  std::vector<int> new_to_old_vertices = construct::linspace(mesh.nv());
  // Utility lambdas to get old, new, and pleat for each vertex.
  auto old = [&](int v) { return new_to_old_vertices[v]; };
  auto new_v = [&](int v) -> std::set<int> & {
    return old_to_new_vertices[old(v)];
  };

  // Returns the index of the origin vertex of the given edge.
  auto fv = [&](auto e) -> int & { return faces[e->fi][e->fvi]; };

  while (!edges_q.empty()) {
    auto edge = edges_q.front();
    edges_q.pop();
    if (!edge || covered_edges.insert(edge).second == false ||
        edge->is_boundary())
      continue;
    // Split the two triangles at the vertex.
    int &v0 = fv(edge);
    // int &v1 = fv(edge->twin()->next());
    if (v0 == fv(edge->twin()->next()) || !has_real_twin(edge->prev())) {
      // Replace v0 with a newly inserted vertex.
      verts.push_back(verts[v0]);
      // Update mappings.
      new_v(v0).insert(verts.size() - 1);
      new_to_old_vertices.push_back(old(v0));
      v0 = verts.size() - 1;

      // Update the adjacent triangle to share v0.
      if (has_real_twin(edge->prev()))
        fv(edge->prev()->twin()) = v0;
    }
    // Split the adjacent edges at the other vertex.
    edges_q.push(edge->prev());
    edges_q.push(edge->next());
    if (!edge->is_boundary()) {
      edges_q.push(edge->twin()->next()->twin());
      edges_q.push(edge->twin()->prev()->twin());
    }
  }

  int num_underlay = verts.size();
  // Add pleat faces.
  // The corresponding pleat index for each original vertex.
  std::vector<int> new_pleats(mesh.nv(), -1);
  double avg_edge_len = mesh.avg_edge_len();
  auto pleat = [&](int v) {
    int &p = new_pleats[old(v)];
    if (p == -1) {
      p = verts.size();
      verts.push_back(verts[old(v)] +
                      mesh.verts[old(v)].normal() * avg_edge_len * 0.5);
    }
    return p;
  };

  for (auto edge : covered_edges) {
    int v0 = fv(edge);
    faces.push_back({v0, pleat(v0), fv(edge->next())});
    if (edge->is_boundary())
      continue;
    int tv0 = fv(edge->twin()->next());
    faces.push_back({pleat(tv0), fv(edge->twin()->next()), fv(edge->twin())});
  }
  for (auto &edge : mesh.edges) {
    if (!edge.is_boundary() || covered_edges.count(&edge) > 0)
      continue;
    faces.push_back({fv(&edge), pleat(fv(edge.next())), fv(edge.next())});
  }

  return fill_pullback_result(orig_mesh, mesh, verts, faces, num_underlay,
                              old_to_new_vertices, new_to_old_vertices,
                              new_pleats, covered_edges);
}

} // namespace cut