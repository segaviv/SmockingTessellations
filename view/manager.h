#pragma once

#include <Eigen/Eigen>
#include <parameterization/param.h>
#include <utils/Hmesh.h>

namespace view {
// Add a triangle mesh.
void add_surface_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      std::string name = "mesh");
// Add a polygon mesh.
void add_surface_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const Eigen::VectorXi &D, std::string name = "remeshed");

void add_surface_mesh(const utils::Hmesh& mesh, std::string name = "mesh");

void add_face_vector_field(std::string mesh_name,
                           const Eigen::MatrixXd &vectors,
                           std::string name = "field");

void show_vector_field(const param::FieldData &field_data,
                       std::string mesh_name = "mesh");

void update_vector_field(const param::FieldData& field_data,
                         std::string mesh_name = "mesh");

void hide_all();
void remove_all();

void show_only(const std::string &name);

void update_visualization();

} // namespace view