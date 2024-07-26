#include "convert.h"
#include "directional/CartesianField.h"
#include <directional/polyvector_to_raw.h>

namespace param {

void multiply_polynomials(Eigen::MatrixXcd &p1, const Eigen::MatrixXcd &p2) {
  assert(p1.rows() == p2.rows());
  Eigen::MatrixXcd newp =
      Eigen::MatrixXcd::Zero(p1.rows(), p1.cols() + p2.cols());
  for (int i = 0; i < p1.cols(); i++)
    for (int j = 0; j < p2.cols(); j++)
      newp.col(i + j).array() += p1.col(i).array() * p2.col(j).array();
  for (int i = 0; i < p1.cols(); i++) {
    newp.col(i + p2.cols()) += p1.col(i);
  }
  for (int i = 0; i < p2.cols(); i++) {
    newp.col(i + p1.cols()) += p2.col(i);
  }

  p1 = newp;
}

Eigen::MatrixXcd raw_to_polyvector(const Eigen::MatrixXd &intField, int N,
                                   bool signSymmetry) {

  Eigen::MatrixXcd pvField =
      Eigen::MatrixXcd::Zero(intField.rows(), intField.cols() / 2);
  Eigen::MatrixXcd actualRoots;
  int actualN;
  if ((signSymmetry) && (N % 2 == 0)) {
    actualRoots.resize(intField.rows(), intField.cols() / 4);
    actualN = N / 2;
    for (int i = 0; i < actualN; i++)
      for (int j = 0; j < intField.rows(); j++)
        actualRoots(j, i) =
            std::complex<double>(intField(j, 2 * i), intField(j, 2 * i + 1));

    // actualRoots=actualRoots*actualRoots;
    actualRoots = actualRoots.array().square();
    // actualRoots = actualRoots.cwiseProduct(actualRoots);
  } else {
    actualRoots.resize(intField.rows(), intField.cols() / 2);
    actualN = N;
    for (int i = 0; i < actualN; i++)
      for (int j = 0; j < intField.rows(); j++)
        actualRoots(j, i) =
            std::complex<double>(intField(j, 2 * i), intField(j, 2 * i + 1));
  }

  int jump = ((signSymmetry) && (N % 2 == 0) ? 2 : 1);
  Eigen::MatrixXcd actualPVField(pvField.rows(), 1);
  actualPVField.col(0) = -actualRoots.col(0);
  for (int i = 1; i < actualN; i++) {
    multiply_polynomials(actualPVField, -actualRoots.col(i));
  }

  for (int i = 0; i < N; i += jump)
    pvField.col(i) = actualPVField.col(i / jump);

  return pvField;
}

Eigen::MatrixXcd
raw_to_polyvector(const directional::CartesianField &raw_field) {
  directional::CartesianField res;
  res.init(*raw_field.tb, directional::fieldTypeEnum::POLYVECTOR_FIELD,
           raw_field.N);
  return raw_to_polyvector(raw_field.intField, raw_field.N, raw_field.N % 2 == 0);
}

directional::CartesianField
polyvector_to_raw(const directional::CartesianField &pv_field,
                  bool signSymmetry) {
  directional::CartesianField res;
  directional::polyvector_to_raw(pv_field, res, signSymmetry);
  return res;
}

} // namespace param