#pragma once

#include "Hmesh.h"

namespace utils {

class ModifiedMesh {
public:
  ModifiedMesh(const utils::Hmesh &mesh) : mesh(mesh) {}

  ModifiedMesh(const utils::Hmesh &mesh, const std::vector<int> &new_to_old)
      : mesh(mesh), new_to_old(new_to_old) {}
  ModifiedMesh(utils::Hmesh &&mesh, std::vector<int> &&new_to_old)
      : mesh(std::move(mesh)), new_to_old(std::move(new_to_old)) {}

  ModifiedMesh &set_original(ModifiedMesh *original) {
    this->original = original;
    return *this;
  }

  // Vertices map.
  std::vector<int> new_to_old;
  utils::Hmesh mesh;

  ModifiedMesh *original = nullptr;
};
} // namespace utils