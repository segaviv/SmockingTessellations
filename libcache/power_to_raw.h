#pragma once

#include <Eigen/Eigen>
#include <directional/CartesianField.h>

namespace directional_ {
void power_to_raw(const directional::CartesianField &powerField, int N,
                  directional::CartesianField &rawField, bool normalize = false);
}