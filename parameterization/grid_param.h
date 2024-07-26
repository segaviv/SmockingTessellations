#pragma once

#include "../utils/Hmesh.h"
#include <Eigen/Eigen>

namespace param {

/**
 * Assuming that the mesh is a pure triangle/quad/hex mesh cut through
 * the singularities, the function returns the initial grid that was pulled back
 * to generate the mesh.
 */
Eigen::MatrixXd grid_param(utils::Hmesh &mesh);

} // namespace param
