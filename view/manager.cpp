#include "manager.h"
#include "glm/fwd.hpp"
#include <parameterization/remeshing.h>
#include "state/state.h"
#include "utils/Hmesh.h"

#include <parameterization/param.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <utils/construct.h>
#include <layout/optimize.h>
#include <utils/conversions.h>
#include <view/smocking.h>

namespace view {
void add_surface_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      std::string name) {
  polyscope::registerSurfaceMesh(name, V, F)->setEnabled(true);
}

void add_surface_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const Eigen::VectorXi &D, std::string name) {
  add_surface_mesh(utils::Hmesh(V, F, D), name);
}

void add_surface_mesh(const utils::Hmesh &mesh, std::string name) {
  auto meshy = polyscope::registerSurfaceMesh(name, mesh.V, mesh.F);
  meshy->setEdgeWidth(1.0);
  meshy->setEnabled(true);
}

void add_face_vector_field(std::string mesh_name,
                           const Eigen::MatrixXd &vectors, std::string name) {
  if (polyscope::getSurfaceMesh(mesh_name) == nullptr)
    return;
  polyscope::getSurfaceMesh(mesh_name)->addFaceVectorQuantity(name, vectors);
}

void show_vector_field(const param::FieldData &field_data,
                       std::string mesh_name) {
  auto mesh = polyscope::getSurfaceMesh(mesh_name);
  if (mesh == nullptr)
    return;
  show_only(mesh_name);
  // Remove all previous fields if they exist.
  for (int i = 0; i < 6; i++) {
    mesh->removeQuantity("field" + std::to_string(i));
  }
  // Add the new field (s).
  for (int i = 0; i < field_data.field.N; i++) {
    auto vf = mesh->addFaceVectorQuantity(
            "field" + std::to_string(i),
            field_data.field.extField.block(
                0, i * 3, field_data.field.extField.rows(), 3));
    vf->setEnabled(true);
    vf->setVectorRadius(0.001, false);
    vf->setVectorLengthScale(0.014, false);
  }
}

void update_vector_field(const param::FieldData &field_data,
                         std::string mesh_name) {
  auto mesh = polyscope::getSurfaceMesh(mesh_name);
  if (mesh == nullptr)
    return;
  for (int i = 0; i < field_data.field.N; i++) {
    auto quantity = (polyscope::SurfaceFaceVectorQuantity *)mesh->getQuantity(
        "field" + std::to_string(i));
    if (!quantity)
      continue;
    quantity->updateData(field_data.field.extField.block(
                0, i * 3, field_data.field.extField.rows(), 3));
    quantity->setEnabled(true);
  }
}

void hide_all() {
  for (auto &structure_type : polyscope::state::structures) {
    for (auto &structure : structure_type.second) {
      structure.second->setEnabled(false);
    }
  }
}
void remove_all() {
  polyscope::removeAllGroups();
  polyscope::removeAllStructures();
}

void show_only(const std::string &name) {
  for (auto &structure_type : polyscope::state::structures) {
    for (auto &structure : structure_type.second) {
      if (structure.first == name)
        structure.second->setEnabled(true);
      else
        structure.second->setEnabled(false);
    }
  }
}

void update_visualization() {
  if (param::optimization_updated) {
    param::optimization_updated = false;
    if (state::field_data.field.N == -2) {
      param::setup_integration(state::field_data.field);
      param::integrate(state::field_data.field);
      auto bla = polyscope::getSurfaceMesh("mesh")->addVertexParameterizationQuantity("param",
      state::field_data.NFunction.leftCols(2) );
      bla->setEnabled(true);
      bla->setCheckerSize(1);
    } else {
      update_vector_field(state::field_data);
    }
  }

  if (layout::layout_updated) {
    layout::layout_updated = false;
    view::show_optimized_layout();
  }
}

} // namespace view