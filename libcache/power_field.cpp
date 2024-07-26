#include "power_field.h"
#include <directional/power_field.h>

namespace directional_ {

void power_field(const directional::TangentBundle &tb,
                 const Eigen::VectorXi &constSpaces,
                 const Eigen::MatrixXd &constVectors,
                 const Eigen::VectorXd &alignWeights, const int N,
                 directional::CartesianField &field) {
  directional::power_field(tb, constSpaces, constVectors, alignWeights, N,
                           field);
}
} // namespace directional_