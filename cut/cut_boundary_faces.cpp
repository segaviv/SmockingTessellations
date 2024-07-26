#include "cut_boundary_faces.h"
#include "utils/Hmesh.h"
#include <set>
#include <utils/find_duplicates.h>

namespace cut {

std::set<int> get_boundary_vertices(const utils::Hmesh& mesh,
                                    const std::set<int>& dup_verts) {
  std::set<int> boundary_vertices;
  for (const auto& edge : mesh.edges) {
    if (edge.is_boundary() && !dup_verts.count(edge.vi) &&
        !dup_verts.count(edge.next()->vi)) {
      boundary_vertices.insert(edge.vi);
      boundary_vertices.insert(edge.next()->vi);
    }
  }
  return boundary_vertices;
}

  utils::Hmesh boundary_faces(const utils::Hmesh &mesh) {
  std::map<int, int> old_to_new_vertex;
  std::vector<int> keep_faces;
  std::set<int> dup_verts = utils::get_duplicated_verts_set(mesh);
  std::set<int> boundary_vertices = get_boundary_vertices(mesh, dup_verts);
  // std::vector<bool> fake_cuts =
  //     get_fake_cut_faces(mesh, v, f, boundary_vertices);

  auto is_boundary_face = [&](int i) {
    for (int v = 0; v < 3; v++) {
      if (boundary_vertices.count(mesh.F[i][v])) {
        return true;
      }
    }
    return false;
  };
  for (int i = 0; i < mesh.faces.size(); i++) {
    // if (!mesh.faces[i].is_boundary() || fake_cuts[i]) {
    if (!is_boundary_face(i)) {
      keep_faces.push_back(i);
      for (int v = 0; v < 3; v++) {
        if (old_to_new_vertex.count(mesh.F[i][v]) == 0) {
          int size = old_to_new_vertex.size();
          old_to_new_vertex[mesh.F[i][v]] = size;
        }
      }
    }
  }
  Eigen::MatrixXd new_verts(old_to_new_vertex.size(), 3);
  for (auto& [old, new_] : old_to_new_vertex) {
    new_verts.row(new_) = mesh.V.row(old);
  }
  Eigen::MatrixXi new_faces(keep_faces.size(), 3);
  for (int i = 0; i < keep_faces.size(); i++) {
    for (int v = 0; v < 3; v++) {
      new_faces(i, v) = old_to_new_vertex[mesh.F[keep_faces[i]][v]];
    }
  }
  return utils::Hmesh(new_verts, new_faces);
}
}