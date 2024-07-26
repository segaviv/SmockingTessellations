#include <Eigen/Eigen>
#include <directional/CartesianField.h>

namespace param {

Eigen::MatrixXcd raw_to_polyvector(const directional::CartesianField &raw_field);

directional::CartesianField polyvector_to_raw(const directional::CartesianField& pv_field,
                                      bool signSymmetry = true);

} // namespace param