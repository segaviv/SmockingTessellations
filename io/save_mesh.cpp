#include "save_mesh.h"
#include <fstream>

namespace io {
void save_raw_mesh_with_uv(const std::string& file_name, const utils::Hmesh& mesh,
                           const Eigen::MatrixXd& uv) {
  std::ofstream output(file_name);
  for (int i = 0; i < mesh.V.rows(); i++) {
    output << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " "
           << (mesh.V.cols() == 3 ? mesh.V(i, 2) : 0.0) << std::endl;
  }
  for (int i = 0; i < mesh.V.rows(); i++) {
    output << "vt " << uv(i, 0) << " " << uv(i, 1) << " " << std::endl;
  }
  for (int i = 0; i < mesh.F.size(); i++) {
    output << "f";
    for (int j = 0; j < mesh.F[i].size(); j++) {
      output << " " << mesh.F[i][j] + 1 << "/" << mesh.F[i][j] + 1;
    }
    output << std::endl;
  }
  output.close();
}

void save_raw_mesh(const std::string& file_name, const utils::Hmesh& mesh) {
  std::ofstream output(file_name);
  for (int i = 0; i < mesh.V.rows(); i++) {
    output << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " "
           << (mesh.V.cols() == 3 ? mesh.V(i, 2) : 0.0) << std::endl;
  }
  for (int i = 0; i < mesh.F.size(); i++) {
    output << "f";
    for (int j = 0; j < mesh.F[i].size(); j++) {
      output << " " << mesh.F[i][j] + 1;
    }
    output << std::endl;
  }
  output.close();
}
}