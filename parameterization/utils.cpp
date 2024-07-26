#include "utils.h"
#include <directional/CartesianField.h>
#include <Optiz/Common/SparseUtils.h>

namespace param {
Eigen::VectorXcd closest_aligned_polyvector(const Eigen::VectorXcd& coefs,
                                            std::complex<double> align_vec,
                                            bool is_even) {
  if (is_even) align_vec = std::pow(align_vec, 2);
  int jump = is_even ? 2 : 1;
  int nv = coefs.size() - 1;
  auto proj_mat = Optiz::sparse<std::complex<double>>(
      {nv}, coefs.size(), [&](int i, auto& x) {
        if (i == 0) return -align_vec * x(0);
        if (i == coefs.size() - 1) return x(coefs.size() - 2);
        return x(i - 1) - align_vec * x(i);
      });
  Eigen::VectorXcd b = Eigen::VectorXcd::Zero(coefs.size());
  b(b.size() - 1) = -align_vec;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver(
      proj_mat.adjoint() * proj_mat);

  return proj_mat * solver.solve(proj_mat.adjoint() * (coefs - b)) + b;
}


void normalize_field(directional::CartesianField& field) {
  Eigen::MatrixXd res = field.intField;
  for (int i = 0; i < res.rows(); i++) {
    for (int j = 0; j < field.N; j++) {
      res.row(i).segment<2>(j * 2) /= res.row(i).segment<2>(j * 2).norm();
    }
  }
  field.set_intrinsic_field(res);
}

}