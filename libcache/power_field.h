#pragma once

#include <directional/CartesianField.h>
#include <Eigen/Eigen>
#include <directional/TangentBundle.h>

namespace directional_ {

void power_field(const directional::TangentBundle& tb,
                                const Eigen::VectorXi& constSpaces,
                                const Eigen::MatrixXd& constVectors,
                                const Eigen::VectorXd& alignWeights,
                                const int N,
                                directional::CartesianField& field);

}