#include "power_to_raw.h"
#include <directional/power_to_raw.h>

namespace directional_ {
void power_to_raw(const directional::CartesianField &powerField, int N,
                  directional::CartesianField &rawField, bool normalize) {
  directional::power_to_raw(powerField, N, rawField, normalize);
}
} // namespace directional_