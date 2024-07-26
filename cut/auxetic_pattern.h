#pragma once

#include "state/state.h"
#include <Eigen/Eigen>
#include <tuple>
#include <utils/Hmesh.h>
#include <parameterization/pullback.h>

namespace cut {
  extern double opening_ang;


param::PullbackResult cut_auxetic_pattern(const utils::Hmesh &mesh);

}