#pragma once
#include <string>
#include <Eigen/Eigen>
#include <utils/Hmesh.h>

namespace io {
void save_raw_mesh_with_uv(const std::string& file_name, const utils::Hmesh& mesh,
                           const Eigen::MatrixXd& uv);

void save_raw_mesh(const std::string& file_name, const utils::Hmesh& mesh);
}