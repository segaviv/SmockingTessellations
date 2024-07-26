#include "find_duplicates.h"

namespace utils { 
  std::map<int, std::set<int>> get_duplicated_verts(const utils::Hmesh& mesh) {
  std::map<int, std::set<int>> duplicated_verts;
  for (int i = 0; i < mesh.verts.size(); i++) {
    for (int j = i + 1; j < mesh.verts.size(); j++) {
      if ((mesh.V.row(i) - mesh.V.row(j)).norm() < 1e-6) {
        duplicated_verts[i].insert(j);
        duplicated_verts[j].insert(i);
      }
    }
  }
  return duplicated_verts;
}

std::set<int> get_duplicated_verts_set(const utils::Hmesh& mesh) {
  std::set<int> duplicated_verts;
  for (int i = 0; i < mesh.verts.size(); i++) {
    for (int j = i + 1; j < mesh.verts.size(); j++) {
      if ((mesh.V.row(i) - mesh.V.row(j)).norm() < 1e-6) {
        duplicated_verts.insert(i);
        duplicated_verts.insert(j);
      }
    }
  }
  return duplicated_verts;
}
}