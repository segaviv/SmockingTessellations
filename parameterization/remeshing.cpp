#include "remeshing.h"
#include <Optiz/Linear/QuadraticObjectiveD.h>
#include <libcache/remeshing.h>
#include <memory>
#include <state/state.h>

namespace param {

std::unique_ptr<directional::IntegrationData> integration_data;
directional::TriMesh meshCut;
directional::CartesianField combedField;

Eigen::VectorXi DPolyMesh;
Eigen::MatrixXi FPolyMesh;
Eigen::MatrixXd VPolyMesh;

double length_ratio = 0.05;

// Integrates the gradient field to get parametrization.
std::shared_ptr<Optiz::QuadraticObjectiveD> gradient_field_integrator;

static Eigen::SparseMatrix<double> grad(const Eigen::MatrixXd &verts,
                                        const Eigen::MatrixXi &faces) {
  const int nf = faces.rows();
  // Iterate over the faces and calculate the coefficients for each
  // face gradient.
  std::vector<Eigen::Triplet<float>> triplets;
  for (int i = 0; i < faces.rows(); i++) {
    const int v1 = faces(i, 0), v2 = faces(i, 1), v3 = faces(i, 2);
    Eigen::Vector3d e1 = verts.row(v2) - verts.row(v1),
                    e2 = verts.row(v3) - verts.row(v2);
    const float e1_len = e1.norm();
    // Coefficients for the first vector of the local frame.
    triplets.push_back(Eigen::Triplet<float>(i, v2, 1. / e1_len));
    triplets.push_back(Eigen::Triplet<float>(i, v1, -1. / e1_len));

    // Coefficients for the second vector of the local frame.
    const Eigen::Vector3d e1_rot = e1.cross(e2).cross(e1).normalized();
    const float coef = e2.dot(e1) / (e1_len * e1_len);
    const float e2_proj_ln = e2.dot(e1_rot);

    triplets.push_back(Eigen::Triplet<float>(i + nf, v3, 1. / e2_proj_ln));
    triplets.push_back(
        Eigen::Triplet<float>(i + nf, v2, -(1. + coef) / e2_proj_ln));
    triplets.push_back(Eigen::Triplet<float>(i + nf, v1, coef / e2_proj_ln));
  }
  Eigen::SparseMatrix<double> grad;
  grad.resize(nf * 2, verts.rows());
  grad.setFromTriplets(triplets.begin(), triplets.end());
  return grad;
}

void setup_integration(const directional::CartesianField &raw_field) {
  directional::CartesianField integrated_field;

  // if (raw_field.N == 2) {
  // gradient_field_integrator = std::make_shared<Optiz::QuadraticObjectiveD>(
  //     std::vector<int>{(int)state::V.rows(), 1}, 2);
  //   gradient_field_integrator->add_quadratic_energy(grad(state::V,
  //   state::F),)
  // }

  if (raw_field.N == 2) {
    integrated_field.init(*raw_field.tb, directional::fieldTypeEnum::RAW_FIELD,
                          4);
    Eigen::MatrixXd field(raw_field.intField.rows(), 8);
    field << raw_field.intField, -raw_field.intField;
    integrated_field.set_intrinsic_field(field);
    integrated_field.matching =
        Eigen::VectorXi::Zero(raw_field.matching.size());
    integrated_field.effort = Eigen::VectorXd::Zero(raw_field.effort.size());

  } else {
    integrated_field = raw_field;
  }

  // Create integration data and set symmetry.
  integration_data =
      std::make_unique<directional::IntegrationData>(integrated_field.N);
  if (integrated_field.N % 3 == 0) {
    integration_data->set_triangular_symmetry(integrated_field.N);
  } else {
    integration_data->set_sign_symmetry(integrated_field.N);
  }

  // Setup integration.
  directional_::setup_integration(integrated_field, *integration_data);
}

void integrate(const directional::CartesianField &raw_field) {
  integration_data->verbose = true;
  integration_data->integralSeamless = true;
  integration_data->roundSeams = false;
  integration_data->lengthRatio = length_ratio;
  directional_::integrate(combedField, *integration_data,
                          state::field_data.NFunction,
                          state::field_data.NCornerFunction);
}

void mesh_isolines(const directional::TriMesh &mesh) {
  directional_::setup_mesh_function_isolines(*integration_data);
  directional_::mesh_function_isolines(mesh, VPolyMesh, DPolyMesh, FPolyMesh);
}

} // namespace param