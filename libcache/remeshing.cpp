#include "remeshing.h"
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/integrate.h>
#include <directional/mesh_function_isolines.h>
#include <directional/setup_integration.h>

namespace directional_ {

directional::MeshFunctionIsolinesData mfiData;
directional::CartesianField combedField;
directional::TriMesh meshCut;

void setup_integration(const directional::CartesianField &field,
                       directional::IntegrationData &intData) {
  directional::setup_integration(field, intData, meshCut, combedField);
}

void integrate(const directional::CartesianField &field,
               directional::IntegrationData &intData,
               Eigen::MatrixXd &NFunction, Eigen::MatrixXd &NCornerFunctions) {
  directional::integrate(combedField, intData, meshCut, NFunction,
                         NCornerFunctions);
}

void setup_mesh_function_isolines(const directional::IntegrationData &intData) {
  mfiData = directional::MeshFunctionIsolinesData();
  directional::setup_mesh_function_isolines(meshCut, intData, mfiData);
}

void mesh_function_isolines(const directional::TriMesh &mesh,
                            Eigen::MatrixXd &VPolyMesh,
                            Eigen::VectorXi &DPolyMesh,
                            Eigen::MatrixXi &FPolyMesh) {
  directional::mesh_function_isolines(mesh, mfiData, /* verbose= */ true,
                                      VPolyMesh, DPolyMesh, FPolyMesh);
}

} // namespace directional_