#pragma once
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif
#include <Eigen/Eigen>
#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

namespace utils {

class Hmesh {
public:
  // Mesh elements.
  struct Edge;
  struct Face;
  struct Vertex;

  Face &face(int i);
  Vertex &vertex(int i);
  Edge &edge(int i);

  struct Face {
    int index;   // index of the face.
    int ei;      // index of one of the edges.
    Hmesh *mesh; // pointer to the mesh.

    Edge *edge() const;

    double area() const;
    Eigen::VectorXd center() const;
    Eigen::VectorXd normal() const;
    bool is_boundary();

    // Following the notation from "Discrete differential operators on polygonal
    // meshes".
    int nf() const { return mesh->F[index].size(); }
  };
  struct Vertex {
    int index;   // index of the vertex.
    int ei;      // index of one of the edges.
    Hmesh *mesh; // pointer to the mesh.

    Eigen::VectorXd coords();
    Eigen::VectorXd normal() const;
    Edge *edge() const;
  };

  struct Edge {
    int index;
    int vi;  // index of the origin vertex.
    int fi;  // index of the face.
    int fvi; // index of the face vertex (vi = F[fi][fvi]).
    int ti;  // index of the twin edge.
    int ni;  // index of the next edge.
    int pi;  // index of the previous edge.
    Hmesh *mesh;

    Vertex *origin() const;
    Vertex *target() const;
    Face *face() const;
    Edge *prev() const;
    Edge *next() const;
    Edge *twin() const;

    double cos_theta();
    double sin_theta();
    double cot_theta();
    double theta();
    Eigen::VectorXd vec();
    double length();
    inline bool is_boundary() const { return !twin() || !face(); }
    std::shared_ptr<void> data;
    template <typename T> T &get_data() {
      return *static_cast<T *>(data.get());
    }
  };

  struct VertexStarIterator {
    VertexStarIterator(Edge *e);
    VertexStarIterator(Edge *e, bool is_begin);
    const VertexStarIterator &operator++();
    bool operator==(const VertexStarIterator &other) const;
    bool operator!=(const VertexStarIterator &other) const;
    Edge *operator*() const { return e; }
    Edge *e;
    bool is_begin;
  };

  struct NodeEdges {
    NodeEdges(const Vertex &v) : e(v.edge()) {}
    VertexStarIterator begin() { return VertexStarIterator(e, true); }
    VertexStarIterator end() { return VertexStarIterator(e, false); }

    Edge *e;
  };

public:
  Eigen::MatrixXd V;
  std::vector<std::vector<int>> F;

  // (optionally cached) Vertex normals.
  Eigen::MatrixXd Vn;

  std::vector<Vertex> verts;
  std::vector<Face> faces;
  std::vector<Edge> edges;
  int total_edges;

  std::vector<std::unordered_map<int, int>> verts_edges;
  Hmesh() = default;
  Hmesh(const Hmesh &other);
  Hmesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
  Hmesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
        const Eigen::VectorXi &D);
  Hmesh(const Eigen::MatrixXd &V, const std::vector<std::vector<int>> &F);
  Hmesh &operator=(const Hmesh &other);

  Hmesh &with_external_boundary_edges();

  Hmesh dual();
  Hmesh remove_unused_verts();

  double avg_edge_len();
  double area() const;

  void compute_vertex_normals();

  std::vector<std::set<int>> vertex_neighbors();

  bool is_boundary_vertex(int index);
  int next_boundary_vert(int index);
  int prev_boundary_vert(int index);
  int manhatten_distance_to_boundary(int index);

  inline int nv() const { return verts.size(); }
  inline int nf() const { return faces.size(); }
};

} // namespace utils