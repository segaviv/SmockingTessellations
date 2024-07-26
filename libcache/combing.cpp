#include "combing.h"
#include <directional/combing.h>

namespace directional_ {
void combing(const directional::CartesianField &rawField,
             directional::CartesianField &combedField) {
  directional::combing(rawField, combedField);
}
} // namespace directional_