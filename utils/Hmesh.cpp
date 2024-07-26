#include "Hmesh.h"

#include <assert.h>

#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <vector>

#include "conversions.h"

namespace utils {

Hmesh::Face &Hmesh::face(int i) {
  assert(i < faces.size() && i >= 0 && "Face index out of bounds");
  return faces[i];
}
Hmesh::Vertex &Hmesh::vertex(int i) {
  assert(i < verts.size() && i >= 0 && "Vertex index out of bounds");
  return verts[i];
}
Hmesh::Edge &Hmesh::edge(int i) {
  assert(i < edges.size() && i >= 0 && "Edge index out of bounds");
  return edges[i];
}

Hmesh::Edge *Hmesh::Face::edge() const { return &mesh->edge(ei); }

Eigen::VectorXd Hmesh::Face::center() const {
  Edge *orig_e = edge();
  Edge *e = orig_e;
  Eigen::VectorXd center = Eigen::VectorXd::Zero(mesh->V.cols());
  int count = 0;
  do {
    center += e->origin()->coords();
    count++;
    e = e->next();
  } while (e != orig_e);
  return center / count;
}

double Hmesh::Face::area() const {
  auto &face = mesh->F[index];
  Eigen::Vector3d vec_area = Eigen::Vector3d::Zero();
  for (int i = 0; i < nf(); i++) {
    vec_area += ((Eigen::Vector3d)mesh->V.row(face[i]))
                    .cross((Eigen::Vector3d)mesh->V.row(face[(i + 1) % nf()]));
  }
  return vec_area.norm() / 2;
}

Eigen::VectorXd Hmesh::Face::normal() const {
  auto &face = mesh->F[index];
  Eigen::Vector3d vec_area = Eigen::Vector3d::Zero();
  for (int i = 0; i < nf(); i++) {
    vec_area += ((Eigen::Vector3d)mesh->V.row(face[i]))
                    .cross((Eigen::Vector3d)mesh->V.row(face[(i + 1) % nf()]));
  }
  return vec_area.normalized();
}

Hmesh::Edge *Hmesh::Vertex::edge() const { return &mesh->edge(ei); }

Eigen::VectorXd Hmesh::Vertex::coords() { return mesh->V.row(index); }
Eigen::VectorXd Hmesh::Vertex::normal() const {
  Eigen::VectorXd normal = Eigen::VectorXd::Zero(mesh->V.cols());
  double total_area = 0;
  Edge *e = edge();
  do {
    normal += e->face()->normal() * e->face()->area();
    total_area += e->face()->area();
    e = e->prev()->twin();
  } while (e != edge() && e != nullptr);
  return (normal / total_area).normalized();
}

bool Hmesh::Face::is_boundary() {
  Edge *orig_e = edge();
  Edge *e = orig_e;
  do {
    if (!e->twin())
      return true;
    e = e->next();
  } while (e != orig_e);
  return false;
}

static std::vector<std::vector<int>> convert_faces(const Eigen::MatrixXi &F,
                                                   const Eigen::VectorXi &D) {
  std::vector<std::vector<int>> faces(F.rows());
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < D(i); j++) {
      faces[i].push_back(F(i, j));
    }
  }
  return faces;
}

Hmesh::Hmesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
             const Eigen::VectorXi &D)
    : Hmesh(V, convert_faces(F, D)) {}

Hmesh::Hmesh(const Eigen::MatrixXd &V, const std::vector<std::vector<int>> &F)
    : V(V), F(F) {
  int nv = V.rows(), nf = F.size();
  verts.reserve(nv);
  faces.reserve(nf);
  int num_edges = 0;
  for (int i = 0; i < F.size(); i++) {
    num_edges += F[i].size();
  }
  edges.reserve(num_edges);
  verts_edges.resize(nv);

  // Push vertices.
  for (int i = 0; i < nv; i++) {
    verts.push_back(Vertex{.index = i, .ei = -1, .mesh = this});
  }

  // Faces and edges.
  int cur_edge_index = 0;
  for (int i = 0; i < nf; i++) {
    // Push face.
    faces.push_back(Face{.index = i, .ei = cur_edge_index, .mesh = this});
    // Push face edges.
    int verts_per_face = F[i].size();
    for (int j = 0; j < verts_per_face; j++) {
      int vi = F[i][j], vip1 = F[i][(j + 1) % verts_per_face];
      int ei = cur_edge_index + j;
      // Edge ei is from vertex vi to vip1.
      edges.push_back(
          Edge{.index = ei,
               .vi = vi,
               .fi = i,
               .fvi = j,
               .ti = -1,
               .ni = cur_edge_index + (j + 1) % verts_per_face,
               .pi = cur_edge_index + (j - 1 + verts_per_face) % verts_per_face,
               .mesh = this});
      // Update vertex edge.
      if (vertex(vi).ei == -1) {
        verts[vi].ei = ei;
      }
      // Update verts_edges.
      verts_edges[vi][vip1] = ei;
      auto it = verts_edges[vip1].find(vi);
      if (it != verts_edges[vip1].end()) {
        edges[ei].ti = it->second;
        edges[it->second].ti = ei;
      }
    }
    cur_edge_index += verts_per_face;
  }
  // Make the boundary edge the first in each vertex.
  for (int i = 0; i < nv; i++) {
    if (verts[i].ei == -1)
      continue;

    Edge *e = verts[i].edge();
    if (!e || !e->twin()) {
      continue;
    }
    do {
      e = e->twin()->next();
    } while (e != verts[i].edge() && e->twin());
    verts[i].ei = e->index;
  }
}

Hmesh &Hmesh::operator=(const Hmesh &other) {
  V = other.V;
  F = other.F;
  verts = other.verts;
  faces = other.faces;
  edges = other.edges;
  verts_edges = other.verts_edges;
  for (auto &v : verts) {
    v.mesh = this;
  }
  for (auto &f : faces) {
    f.mesh = this;
  }
  for (auto &e : edges) {
    e.mesh = this;
  }
  return *this;
}

inline std::vector<std::vector<int>> to_vector(const Eigen::MatrixXi &F) {
  std::vector<std::vector<int>> F_vec(F.rows());
  for (int i = 0; i < F.rows(); i++) {
    F_vec[i].resize(F.cols());
    for (int j = 0; j < F.cols(); j++) {
      F_vec[i][j] = F(i, j);
    }
  }
  return std::move(F_vec);
}

Hmesh::Hmesh(const Hmesh &other)
    : V(other.V), F(other.F), verts(other.verts), faces(other.faces),
      edges(other.edges), verts_edges(other.verts_edges) {
  for (auto &v : verts) {
    v.mesh = this;
  }
  for (auto &f : faces) {
    f.mesh = this;
  }
  for (auto &e : edges) {
    e.mesh = this;
  }
}

Hmesh::Hmesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    : Hmesh(V, to_vector(F)) {}

Hmesh &Hmesh::with_external_boundary_edges() {
  for (Edge &e : edges) {
    if (!e.twin()) {
      edges.push_back(Edge{.index = (int)edges.size(),
                           .vi = e.next()->vi,
                           .fi = -1,
                           .ti = e.index,
                           .ni = -1,
                           .pi = -1});
      e.ti = edges.size() - 1;
      verts_edges[e.twin()->vi][e.vi] = e.twin()->index;
    }
  }
  return *this;
}

Hmesh Hmesh::dual() {
  Eigen::MatrixXd new_v(F.size(), 3);
  std::vector<std::vector<int>> new_faces;
  for (int i = 0; i < faces.size(); i++) {
    new_v.row(i) = faces[i].center();
  }
  for (auto &vert : verts) {
    std::vector<int> new_face;
    for (auto e : Hmesh::NodeEdges(vert)) {
      new_face.push_back(e->fi);
    }
    if (new_face.size() > 2)
      new_faces.push_back(new_face);
  }
  // for (int i = 0; i < V.rows(); i++) {
  //   new_v.row(i + faces.size()) = V.row(i);
  // }
  // for (int i = 0; i < verts.size(); i++) {
  //   Edge* e = verts[i].edge;
  //   if (!e || !e->twin) {
  //     continue;
  //   }
  //   do {
  //     new_faces.push_back({e->face->index, e->twin->face->index});
  //     e = e->twin->next;
  //   } while (e != verts[i].edge && e->twin);
  // }
  return Hmesh(new_v, new_faces);
}

Hmesh Hmesh::remove_unused_verts() {
  std::map<int, int> new_vert_map;
  std::vector<int> new_to_old;
  std::vector<std::vector<int>> faces = F;
  for (auto &face : faces) {
    for (int j = 0; j < face.size(); j++) {
      auto res = new_vert_map.insert({face[j], new_to_old.size()});
      if (res.second)
        new_to_old.push_back(face[j]);
      face[j] = res.first->second;
    }
  }
  return Hmesh(V(new_to_old, Eigen::all), faces);
}

double Hmesh::avg_edge_len() {
  double sum = 0;
  int count = 0;
  for (Edge &e : edges) {
    if (!e.twin() || e.vi < e.next()->vi) {
      sum += e.vec().norm();
      count++;
    }
  }
  return sum / count;
}

double Hmesh::area() const {
  double sum = 0;
  for (auto &f : faces) {
    sum += f.area();
  }
  return sum;
}

bool Hmesh::is_boundary_vertex(int index) {
  for (auto ei : verts_edges[index]) {
    auto e = edges[ei.second];
    if (e.twin() == nullptr || e.twin()->fi == -1) {
      return true;
    }
  }
  return false;
}

int Hmesh::next_boundary_vert(int index) {
  return verts[index].edge()->next()->vi;
}

int Hmesh::prev_boundary_vert(int index) {
  auto e = verts[index].edge();
  while (!e->prev()->is_boundary())
    e = e->prev()->twin();
  return e->prev()->vi;
}

Hmesh::VertexStarIterator::VertexStarIterator(Edge *e)
    : VertexStarIterator(e, true) {}
Hmesh::VertexStarIterator::VertexStarIterator(Edge *e, bool is_begin)
    : e(e), is_begin(is_begin) {}
const Hmesh::VertexStarIterator &Hmesh::VertexStarIterator::operator++() {
  e = e->prev()->twin();
  is_begin = false;
  return *this;
}
bool Hmesh::VertexStarIterator::operator==(
    const VertexStarIterator &other) const {
  return !e || !e->face() || (e == other.e && is_begin == other.is_begin);
}
bool Hmesh::VertexStarIterator::operator!=(
    const VertexStarIterator &other) const {
  return !(operator==(other));
}

Hmesh::Vertex *Hmesh::Edge::origin() const { return &mesh->vertex(vi); }
Hmesh::Face *Hmesh::Edge::face() const { return &mesh->face(fi); }
Hmesh::Edge *Hmesh::Edge::prev() const { return &mesh->edge(pi); }
Hmesh::Edge *Hmesh::Edge::next() const { return &mesh->edge(ni); }
Hmesh::Edge *Hmesh::Edge::twin() const {
  if (ti == -1)
    return nullptr;
  return &mesh->edge(ti);
}

double Hmesh::Edge::cos_theta() {
  return -next()->vec().normalized().dot(prev()->vec().normalized());
}
double Hmesh::Edge::sin_theta() {
  return static_cast<Eigen::Vector3d>(next()->vec().normalized())
      .cross(static_cast<Eigen::Vector3d>(prev()->vec().normalized()))
      .norm();
}
double Hmesh::Edge::cot_theta() {
  // return cos_theta() / sin_theta();
  Eigen::Vector3d prev_vec = prev()->vec();
  Eigen::Vector3d next_vec = next()->vec();
  return -prev_vec.dot(next_vec) / prev_vec.cross(next_vec).norm();
}
double Hmesh::Edge::theta() { return acos(cos_theta()); }
Hmesh::Vertex *Hmesh::Edge::target() const { return next()->origin(); }
Eigen::VectorXd Hmesh::Edge::vec() {
  return next()->origin()->coords() - origin()->coords();
}
double Hmesh::Edge::length() { return vec().norm(); }

std::vector<std::set<int>> Hmesh::vertex_neighbors() {
  std::vector<std::set<int>> neighbors(V.rows());
  for (int i = 0; i < F.size(); i++) {
    for (int j = 0; j < F[i].size(); j++) {
      int v0 = F[i][j], v1 = F[i][(j + 1) % F[i].size()];
      neighbors[v0].insert(v1);
      neighbors[v1].insert(v0);
    }
  }
  return neighbors;
}

int Hmesh::manhatten_distance_to_boundary(int index) {
  std::queue<std::pair<int, int>> q;
  std::set<int> visited;
  if (is_boundary_vertex(index)) {
    return 0;
  }
  q.push({index, 0});
  while (!q.empty()) {
    auto [cur, dist] = q.front();
    q.pop();
    if (is_boundary_vertex(cur)) {
      return dist;
    }
    if (visited.count(cur)) {
      continue;
    }
    visited.insert(cur);
    for (auto edge : NodeEdges(verts[cur])) {
      if (!visited.count(edge->ti)) {
        q.push({edge->ti, dist + 1});
      }
    }
  }
  return -1;
}

void Hmesh::compute_vertex_normals() {
  Vn = Eigen::MatrixXd::Zero(V.rows(), V.cols());
  for (int i = 0; i < F.size(); i++) {
    auto face_normal = faces[i].normal();
    for (int j = 0; j < F[i].size(); j++) {
      Vn(F[i][j], 0) += face_normal(0);
      Vn(F[i][j], 1) += face_normal(1);
      Vn(F[i][j], 2) += face_normal(2);
    }
  }
  Vn.colwise().normalize();
}

} // namespace utils