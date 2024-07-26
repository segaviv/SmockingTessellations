#include "load_mesh.h"
#include <utils/conversions.h>
#include <fstream>

namespace io {
utils::Hmesh load_mesh(const std::string& filename) {
  std::ifstream input(filename);
  std::string line;
  std::vector<Eigen::Vector3d> verts;
  std::vector<std::vector<int>> faces;
  while (std::getline(input, line)) {
    if (line[0] == 'v' && line[1] == ' ') {
      std::stringstream ss(line);
      std::string v;
      double x, y, z;
      ss >> v >> x >> y >> z;
      verts.push_back(Eigen::Vector3d(x, y, z));
    } else if (line[0] == 'f') {
      std::stringstream ss(line);
      std::string f;
      ss >> f;
      int v;
      std::vector<int> face;
      while (ss >> v) {
        face.push_back(v - 1);
        char bla;
        while (ss >> std::noskipws >> bla && !std::isspace(bla)) {
        }
      }
      faces.push_back(face);
    }
  }
  return utils::Hmesh(convert::to_eig_mat(verts), faces);
}
}