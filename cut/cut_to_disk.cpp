#include "cut_to_disk.h"
#include "utils/ModifiedMesh.h"
#include "utils/conversions.h"

#include <iostream>
#include <map>
#include <queue>
#include <state/state.h>
#include <tuple>
#include <utils/slicing.h>

namespace cut {

std::vector<int> singularity_indices(utils::Hmesh &mesh, int symmetry) {
  assert((symmetry % 3 == 0 || symmetry % 4 == 0) && "Symmetry must be 3 or 4");
  int valency = (symmetry % 3 == 0) ? 6 : 4;
  return slice::indices(mesh.nv(), [&](int i) {
    return !mesh.is_boundary_vertex(i) && mesh.verts_edges[i].size() != valency;
  });
}

std::vector<int> connect_boundaries(utils::Hmesh &mesh,
                                    std::vector<std::set<int>> &boundaries,
                                    int index,
                                    const std::set<int> &avoid_verts) {
  // Set targets.
  std::set<int> targets;
  for (int i = 0; i < boundaries.size(); i++) {
    if (i == index)
      continue;
    for (auto &v : boundaries[i]) {
      if (!avoid_verts.count(v)) {
        targets.insert(v);
      }
    }
  }
  // Set sources.
  std::queue<int> queue;
  // Parent id and 'straightness'.
  std::map<int, std::pair<int, double>> parent;
  for (auto &v : boundaries[index]) {
    // if (!dup_verts.count(v)) {
    queue.push(v);
    // }
    parent[v] = {-1, 0.0};
  }
  auto vertex_neighbors = mesh.vertex_neighbors();
  // BFS until reaching one of the targets and return the path.
  while (!queue.empty()) {
    int current = queue.front();
    queue.pop();
    if (targets.find(current) != targets.end()) {
      // Found a path.
      std::vector<int> path;
      while (current != -1) {
        path.push_back(current);
        current = parent[current].first;
      }
      return path;
    }

    // Find unvisited neighbors.
    std::vector<std::pair<int, double>> unvisited_neighbors;
    int cur_parrent = parent[current].first;
    for (auto &n : vertex_neighbors[current]) {
      if (avoid_verts.count(n))
        continue;
      double straightness =
          cur_parrent == -1
              ? 0.0
              : (mesh.V.row(current) - mesh.V.row(cur_parrent))
                    .normalized()
                    .dot((mesh.V.row(n) - mesh.V.row(current)).normalized());
      auto p = parent.find(n);
      if (p != parent.end())
        continue;
      unvisited_neighbors.push_back({n, straightness});
    }
    std::sort(
        unvisited_neighbors.begin(), unvisited_neighbors.end(),
        [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
          return a.second > b.second;
        });

    for (auto &n : unvisited_neighbors) {
      if (parent.find(n.first) == parent.end()) {
        queue.push(n.first);
      }
      parent[n.first] = {current, n.second};
    }

    // Push neighbors.
    // for (auto& n : vertex_neighbors[current]) {
    //   if (parent.find(n) == parent.end()) {
    //     parent[n] = current;
    //     queue.push(n);
    //   }
    // }
  }
  return std::vector<int>();
}

std::vector<std::set<int>>
boundary_loops(utils::Hmesh &mesh,
               const std::set<int> &additional_boundary_verts) {
  std::set<int> boundary_verts;
  for (int i = 0; i < mesh.verts.size(); i++) {
    if (mesh.is_boundary_vertex(i)) {
      boundary_verts.insert(i);
    }
  }
  std::vector<std::set<int>> boundaries;
  while (!boundary_verts.empty()) {
    std::set<int> boundary;
    int current = *boundary_verts.begin();
    while (boundary.count(current) == 0) {
      boundary.insert(current);
      boundary_verts.erase(current);
      if (mesh.verts[current].edge()->twin()) {
        std::cout << "Not real boundary : " << current << std::endl;
        // Not a real boundary.
        break;
      }
      current = mesh.verts[current].edge()->next()->vi;
    }
    boundaries.push_back(boundary);
  }
  for (int v : additional_boundary_verts) {
    bool add_v = true;
    for (auto &b : boundaries) {
      if (b.count(v)) {
        add_v = false;
        break;
      }
    }
    if (add_v) {
      boundaries.push_back({v});
    }
  }
  return boundaries;
}

static void update_face_vertex_index(std::vector<int> &face, int old_v,
                                     int new_v) {
  for (int i = 0; i < face.size(); i++) {
    if (face[i] == old_v) {
      face[i] = new_v;
      return;
    }
  }
}

static std::tuple<Eigen::MatrixXd, std::vector<std::vector<int>>>
make_cut(utils::Hmesh &mesh, const std::vector<int> &path,
         std::set<int> &avoid_verts, std::vector<int> &new_to_old) {
  int is_start_boundary = mesh.is_boundary_vertex(path[0]) ? 1 : 0;
  int is_end_boundary = mesh.is_boundary_vertex(path.back()) ? 1 : 0;
  int new_verts = is_start_boundary + is_end_boundary + path.size() - 2;
  Eigen::MatrixXd new_V(mesh.V.rows() + new_verts, 3);
  new_V.topRows(mesh.V.rows()) = mesh.V;
  std::vector<std::vector<int>> new_faces = mesh.F;

  int new_V_idx = mesh.V.rows();
  for (int i = 1 - is_start_boundary; i < path.size() - 1 + is_end_boundary;
       i++) {
    // Update faces.
    auto cur_edge = i > 0
                        ? &mesh.edge(mesh.verts_edges[path[i]].at(path[i - 1]))
                        : mesh.verts[path[i]].edge();
    auto end_edge = i < path.size() - 1
                        ? &mesh.edge(mesh.verts_edges[path[i]].at(path[i + 1]))
                        : nullptr;
    while (cur_edge != end_edge) {
      update_face_vertex_index(new_faces[cur_edge->fi], cur_edge->vi,
                               new_V_idx);
      cur_edge = cur_edge->prev()->twin();
    }
    // Update vertices.
    new_V.row(new_V_idx) = mesh.V.row(path[i]);
    new_to_old.push_back(path[i]);
    avoid_verts.insert(new_V_idx);
    new_V_idx++;
  }
  return std::make_tuple(new_V, new_faces);
}

utils::ModifiedMesh to_disk(utils::Hmesh &mesh, int symmetry,
                            bool add_cut_twin_edges) {
  std::cout << "Cutting to disk... ";
  utils::Hmesh result(mesh.V, mesh.F);
  std::vector<int> new_to_old(mesh.V.rows());
  for (int i = 0; i < mesh.verts.size(); i++) {
    new_to_old[i] = i;
  }

  std::set<int> boundary_verts =
      convert::to<std::set<int>>(singularity_indices(mesh, symmetry));
  std::set<int> avoid_verts = boundary_verts;
  auto boundaries = boundary_loops(result, boundary_verts);
  while (boundaries.size() > 1) {
    std::cout << "Number of connected components: " << boundaries.size()
              << std::endl;
    std::vector<int> path;
    // Try to find a path to connect boundary 0 with the rest.
    avoid_verts.clear();
    path = connect_boundaries(result, boundaries, boundaries.size() - 1,
                              avoid_verts);
    if (path.size() == 0) {
      std::cout << "Failed to cut to disk. Is there more than one connected "
                   "components?"
                << std::endl;
      return utils::ModifiedMesh(result, new_to_old);
    }
    avoid_verts.insert(path.begin(), path.end());
    std::cout << "Making cut... ";
    auto [new_V, new_faces] = make_cut(result, path, avoid_verts, new_to_old);
    result = utils::Hmesh(new_V, new_faces);
    boundaries = boundary_loops(result, boundary_verts);
  }
  std::cout << "Number of connected components: " << boundaries.size()
            << std::endl;
  return utils::ModifiedMesh(result, new_to_old);
}
} // namespace cut
