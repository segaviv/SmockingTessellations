

#include <Eigen/Eigen>
#include <directional/CartesianField.h>

namespace param {
Eigen::VectorXcd closest_aligned_polyvector(const Eigen::VectorXcd &coefs,
                                            std::complex<double> align_vec,
                                            bool is_even);

void normalize_field(directional::CartesianField &field);
} // namespace param