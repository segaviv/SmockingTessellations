#include "state.h"
#include "igl/doublearea.h"
#include <igl/face_areas.h>
#include <igl/read_triangle_mesh.h>
#include <io/load_mesh.h>

namespace state {
utils::Hmesh mesh, remeshed, remeshed_cut;
int symmetry = 0;
param::FieldData field_data;
param::PullbackResult pullback_result;
Eigen::MatrixXd optimized_layout;
directional::CartesianField initial_field;

std::unique_ptr<visualization::SeamlessARAP> pleats_arap, underlay_arap;
Eigen::MatrixXd arap_verts;
Eigen::MatrixXi arap_faces;

int rep_x = 1, rep_y = 1;
float stitching_scale = 1.0;
smocking::Tangram tangram =
    smocking::Tangram(smocking::UnitSmockingPattern::ARROW);

void reset() {
  field_data = param::FieldData();
  initial_field = directional::CartesianField();
  pullback_result = param::PullbackResult();
  optimized_layout = Eigen::MatrixXd();
  state::remeshed = utils::Hmesh();
  state::remeshed_cut = utils::Hmesh();

  pleats_arap.reset();
  underlay_arap.reset();
  arap_verts = Eigen::MatrixXd();
  arap_faces = Eigen::MatrixXi();
}

bool load_mesh(const std::string &filename) {
  state::mesh = io::load_mesh(filename);
  double area = state::mesh.area();
  mesh.V /= sqrt(area);
  mesh.V.rowwise() -= mesh.V.colwise().mean();
  mesh.compute_vertex_normals();
  reset();
  return true;
}

void init_vector_field(int symmetry, bool align_to_boundary) {
  Eigen::MatrixXi F(state::mesh.F.size(), 3);
  for (int i = 0; i < state::mesh.F.size(); i++) {
    for (int j = 0; j < 3; j++) {
      F(i, j) = state::mesh.F[i][j];
    }
  }
  param::init_field(state::mesh.V, F, field_data, symmetry, align_to_boundary);
  initial_field = field_data.field;
}

void set_unit_pattern(smocking::UnitSmockingPattern pattern) {
  tangram = smocking::Tangram(pattern);
}

} // namespace state