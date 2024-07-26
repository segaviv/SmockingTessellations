#pragma once

#include <directional/CartesianField.h>
#include <directional/setup_integration.h>

namespace directional_ {
void setup_integration(const directional::CartesianField &field,
                       directional::IntegrationData &intData);

void integrate(const directional::CartesianField &field,
               directional::IntegrationData &intData,
               Eigen::MatrixXd &NFunction, Eigen::MatrixXd &NCornerFunctions);

void setup_mesh_function_isolines(const directional::TriMesh &meshCut,
                                  const directional::IntegrationData &intData);

void setup_mesh_function_isolines(const directional::IntegrationData &intData);

void mesh_function_isolines(const directional::TriMesh &mesh,
                            Eigen::MatrixXd &VPolyMesh,
                            Eigen::VectorXi &DPolyMesh,
                            Eigen::MatrixXi &FPolyMesh);
} // namespace directional_