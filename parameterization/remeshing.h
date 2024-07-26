#pragma once

#include <directional/CartesianField.h>

namespace param {

extern Eigen::VectorXi DPolyMesh;
extern Eigen::MatrixXi FPolyMesh;
extern Eigen::MatrixXd VPolyMesh;

extern double length_ratio;

void setup_integration(const directional::CartesianField &raw_field);
void integrate(const directional::CartesianField &raw_field);

void mesh_isolines(const directional::TriMesh &mesh);
} // namespace param