#include "grid_param.h"

#include <iostream>
#include <queue>
#include <set>

#include <utils/Hmesh.h>

namespace param {

static void parameterize_face_by_edge(utils::Hmesh& mesh,
                                      utils::Hmesh::Edge* edge,
                                      double avg_edge_len, int degree,
                                      Eigen::MatrixXd& UV) {
  double angle = (degree - 2) * M_PI / degree;
  Eigen::Rotation2Dd rot(M_PI - angle);

  Eigen::Vector2d uv_edge =
      (UV.row(edge->next()->vi) - UV.row(edge->vi)).normalized() * avg_edge_len;
  auto cur_edge = edge->next();
  while (cur_edge->next() != edge) {
    uv_edge = rot * uv_edge;
    Eigen::Vector2d new_uv =
        ((Eigen::Vector2d)UV.row(cur_edge->vi)) + uv_edge;
    UV.row(cur_edge->next()->vi) = new_uv;
    cur_edge = cur_edge->next();
  }
}

static int get_mesh_common_degree(utils::Hmesh& mesh) {
  std::map<int, int> historgram;
  for (int i = 0; i < mesh.F.size(); i++) {
    historgram[mesh.F[i].size()]++;
  }
  int max_count = 0;
  int max_degree = 0;
  for (auto& [degree, count] : historgram) {
    if (count > max_count) {
      max_count = count;
      max_degree = degree;
    }
  }
  return max_degree;
}

Eigen::MatrixXd grid_param(utils::Hmesh& mesh) {
  double avg_edge_len = mesh.avg_edge_len();
  avg_edge_len = 1;

  Eigen::MatrixXd UV = Eigen::MatrixXd::Zero(mesh.V.rows(), 2);
  std::queue<utils::Hmesh::Edge*> Q;
  std::set<int> visited_faces;
  while (visited_faces.size() < mesh.F.size()) {
    UV.row(mesh.faces[0].edge()->next()->vi) = Eigen::Vector2d(avg_edge_len, 0);
    Q.push(mesh.faces[0].edge());
    while (!Q.empty()) {
      utils::Hmesh::Edge* edge = Q.front();
      Q.pop();
      if (visited_faces.find(edge->fi) != visited_faces.end()) {
        continue;
      }
      visited_faces.insert(edge->fi);
      int f_degree = mesh.F[edge->fi].size();
      parameterize_face_by_edge(mesh, edge, avg_edge_len, f_degree, UV);
      utils::Hmesh::Edge* cur_edge = edge;
      do {
        if (cur_edge->twin()) Q.push(cur_edge->twin());
        cur_edge = cur_edge->next();
      } while (cur_edge != edge);
    }
  }

  return UV;
}
}  // namespace param
