#pragma once

#include <map>
#include <set>

#include "Hmesh.h"

namespace utils { 
  std::map<int, std::set<int>> get_duplicated_verts(const utils::Hmesh& mesh);
  std::set<int> get_duplicated_verts_set(const utils::Hmesh& mesh);
}