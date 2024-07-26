
#include "param.h"
#include "Optiz/Linear/LinearExpressionVec2.h"
#include "convert.h"
#include "directional/CartesianField.h"
#include "directional/TangentBundle.h"
#include "directional/TriMesh.h"
#include "libcache/combing.h"
#include "libcache/principal_matching.h"
#include "utils.h"
#include "utils/Hmesh.h"
#include <Optiz/Linear/QuadraticObjectiveCD.h>
#include <Optiz/Linear/QuadraticObjectiveD.h>
#include <Optiz/Problem.h>
#include <chrono>
#include <complex>
#include <libcache/power_field.h>
#include <libcache/power_to_raw.h>
#include <memory>
#include <unordered_map>

namespace param {

// Parameters for optimization.
struct OptimizationParams {
  int iter = 0;
  double smoothness_weight = 10.0;
  double alignment_weight = 100.0;
  double closeness_weight = 1e-1;
  double rosy_weight = 1e-1;
  bool operator==(const OptimizationParams &other) const {
    return smoothness_weight == other.smoothness_weight &&
           alignment_weight == other.alignment_weight &&
           closeness_weight == other.closeness_weight &&
           rosy_weight == other.rosy_weight;
  }
  bool operator!=(const OptimizationParams &other) const {
    return !(*this == other);
  }
};
OptimizationParams params, cached_params{.smoothness_weight = -1.0};
Eigen::VectorXi cached_matching;
std::unique_ptr<Optiz::QuadraticObjectiveD> reduce_curl_objective;
std::unique_ptr<Optiz::QuadraticObjectiveCD> smooth_polyvector_objective;

bool stop_optimization = true;
bool optimization_updated;

// Alignment.
Eigen::VectorXi const_faces;
Eigen::MatrixXd const_vectors;

static void setup_constrainted_faces(const Eigen::MatrixXd &V,
                                     const Eigen::MatrixXi &F,
                                     bool align_to_boundary) {
  if (!align_to_boundary) {
    const_faces = Eigen::VectorXi{0};
    const_vectors = Eigen::MatrixXd(1, 3);
    const_vectors << (V.row(F(0, 1)) - V.row(F(0, 0))).normalized();
    return;
  }

  // Find boundary edges, and add constraints for the adjacent faces.
  utils::Hmesh mesh(V, F);
  std::unordered_map<int, Eigen::VectorXd> constraints;
  for (auto &edge : mesh.edges) {
    if (edge.is_boundary()) {
      constraints[edge.fi] = edge.vec().normalized();
    }
  }
  const_faces = Eigen::VectorXi(constraints.size());
  const_vectors = Eigen::MatrixXd(constraints.size(), 3);
  int i = 0;
  for (auto &[face, vec] : constraints) {
    const_faces[i] = face;
    const_vectors.row(i) = vec;
    i++;
  }
}

static void add_orthogonal_field(directional::CartesianField &field) {
  Eigen::MatrixXd new_int_field(field.intField.rows(), 4);
  new_int_field << field.intField,
      field.intField * Eigen::Rotation2Dd(-M_PI / 2).toRotationMatrix();
  field.N = 2;
  field.set_intrinsic_field(new_int_field);
}

void init_field(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                FieldData &res, int symmetry, bool align_to_boundary) {
  setup_constrainted_faces(V, F, align_to_boundary);

  res.mesh.set_mesh(V, F);
  res.tangent_bundle.init(res.mesh);

  directional::CartesianField power_field;
  directional_::power_field(res.tangent_bundle, const_faces, const_vectors,
                            Eigen::VectorXd::Ones(const_faces.size()), symmetry,
                            power_field);
  directional::CartesianField tmp;
  directional_::power_to_raw(power_field, symmetry, tmp, true);
  directional_::principal_matching(tmp);
  directional_::combing(tmp, res.field);

  // If symmetry is 1, add the field rotated by pi/2 as the second field.
  if (symmetry == 1) {
    add_orthogonal_field(res.field);
  }
  cached_params.alignment_weight = -1;
  params = OptimizationParams();
}

static std::tuple<std::complex<double>, std::complex<double>>
get_shared_edge_local_frames(directional::IntrinsicFaceTangentBundle *tb,
                             int i) {
  const directional::TriMesh *mesh = tb->mesh;
  int f0 = tb->adjSpaces(i, 0), f1 = tb->adjSpaces(i, 1);

  Eigen::Vector3d shared_edge =
      (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
  std::complex<double> ef0(shared_edge.dot(mesh->FBx.row(f0)),
                           shared_edge.dot(mesh->FBy.row(f0)));
  std::complex<double> ef1(shared_edge.dot(mesh->FBx.row(f1)),
                           shared_edge.dot(mesh->FBy.row(f1)));
  return {ef0, ef1};
}

static std::tuple<Eigen::Vector2d, Eigen::Vector2d>
get_shared_edge_local_frames_real(directional::IntrinsicFaceTangentBundle *tb,
                                  int i) {
  const directional::TriMesh *mesh = tb->mesh;
  int f0 = tb->adjSpaces(i, 0), f1 = tb->adjSpaces(i, 1);

  Eigen::Vector3d shared_edge =
      (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
  Eigen::Vector2d ef0(shared_edge.dot(mesh->FBx.row(f0)),
                      shared_edge.dot(mesh->FBy.row(f0)));
  Eigen::Vector2d ef1(shared_edge.dot(mesh->FBx.row(f1)),
                      shared_edge.dot(mesh->FBy.row(f1)));
  return {ef0, ef1};
}

static double get_edge_weight(directional::IntrinsicFaceTangentBundle *tb,
                              int i) {
  int f0 = tb->adjSpaces(i, 0), f1 = tb->adjSpaces(i, 1);
  double f0_area = tb->mesh->faceAreas(f0), f1_area = tb->mesh->faceAreas(f1);
  double edge_len = (tb->mesh->V.row(tb->mesh->EV(i, 1)) -
                     tb->mesh->V.row(tb->mesh->EV(i, 0)))
                        .norm();
  return 3 * (edge_len * edge_len) / (f0_area + f1_area);
}

static std::complex<double>
get_face_constraint(directional::IntrinsicFaceTangentBundle *tb, int i) {
  int f = const_faces(i);
  std::complex<double> rotated_vec(
      const_vectors.row(i).dot(tb->mesh->FBy.row(f)),
      -const_vectors.row(i).dot(tb->mesh->FBx.row(f)));
  return rotated_vec;
}
static Eigen::Vector2d
get_face_constraint_real(directional::IntrinsicFaceTangentBundle *tb, int i) {
  int f = const_faces(i);
  Eigen::Vector2d rotated_vec(const_vectors.row(i).dot(tb->mesh->FBy.row(f)),
                              -const_vectors.row(i).dot(tb->mesh->FBx.row(f)));
  return rotated_vec;
}

static Eigen::VectorXcd
get_closest_aligned_polyvector(directional::IntrinsicFaceTangentBundle *tb,
                               const Eigen::MatrixXcd &pv_field, int i) {
  std::complex<double> u = get_face_constraint(tb, i);
  int f = const_faces(i);
  int N = pv_field.cols();
  bool sign_symmetry = pv_field.cols() % 2 == 0;
  int nc = sign_symmetry ? N / 2 : N;
  Eigen::VectorXcd coefs(nc);
  if (sign_symmetry) {
    for (int i = 0; i < N - 2; i += 2) {
      coefs(i / 2) = pv_field(f, i);
    }
    // coefs = pv_field(f, Eigen::seq(0, N - 2, 2));
  } else {
    coefs = pv_field.row(f);
  }
  return closest_aligned_polyvector(coefs, u, sign_symmetry);
}

static void
compute_smooth_pv_objective(directional::IntrinsicFaceTangentBundle *tb,
                            const Eigen::MatrixXcd &pv_field) {
  int N = pv_field.cols();
  bool sign_symmetry = N % 2 == 0;
  // Number of coefficients per face.
  int nc = sign_symmetry ? N / 2 : N;
  int nf = pv_field.rows();
  int jump = sign_symmetry ? 2 : 1;

  smooth_polyvector_objective =
      std::make_unique<Optiz::QuadraticObjectiveCD>(std::vector<int>{nf, nc});
  auto start = std::chrono::steady_clock::now();
  // Smoothness term.
  smooth_polyvector_objective->add_weighted_equations(
      tb->adjSpaces.rows(), [&](int i, auto &x, const auto &add_equation) {
        int f = tb->adjSpaces(i, 0), g = tb->adjSpaces(i, 1);
        if (f == -1 || g == -1)
          return;
        double weight = params.smoothness_weight * get_edge_weight(tb, i);
        // Get the shared edge in local coords of faces f, g.
        auto [ef, eg] = get_shared_edge_local_frames(tb, i);
        for (int n = 0; n < nc; n++) {
          // e^(N-n).
          auto efn = std::pow(std::conj(ef), N - jump * n);
          auto egn = std::pow(std::conj(eg), N - jump * n);
          // nth coefficient of the polyvector field.
          add_equation(weight, efn * x(f, n) - egn * x(g, n));
        }
      });

  // Closeness term.
  smooth_polyvector_objective->add_weighted_equations_with_rhs(
      nf, [&](int i, auto &x, const auto &add_equation) {
        double face_area = tb->mesh->faceAreas(i);
        double weight = params.closeness_weight * face_area;
        for (int n = 0; n < nc; n++) {
          add_equation(weight, x(i, n), pv_field(i, n));
        }
      });

  // Orthogonality term (try to make it n-rosy).
  smooth_polyvector_objective->add_weighted_equations_with_rhs(
      nf, [&](int i, auto &x, const auto &add_equation) {
        double weight = params.rosy_weight * tb->mesh->faceAreas(i);
        std::complex<double> normalized =
            pv_field(i, 0) / std::abs(pv_field(i, 0));
        // Try to make it as close to z^N = e^(it) as possible.
        add_equation(weight, x(i, 0), normalized);
      });
  smooth_polyvector_objective->add_weighted_equations(
      nf, [&](int i, auto &x, const auto &add_equation) {
        double weight = params.rosy_weight * tb->mesh->faceAreas(i);
        for (int n = 1; n < nc; n++) {
          add_equation(weight, x(i, n));
        }
      });

  // Alignment term.
  if (const_faces.size() > 1) {
    smooth_polyvector_objective->add_weighted_equations_with_rhs(
        const_faces.size(), [&](int i, auto &x, const auto &add_equation) {
          int f = const_faces(i);
          double weight = params.alignment_weight * tb->mesh->faceAreas(f);
          auto closest = get_closest_aligned_polyvector(tb, pv_field, i);
          for (int n = 0; n < nc; n++) {
            add_equation(weight, x(f, n), closest(n));
          }
        });
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "smooth polyvector objective time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
}

static void
update_smooth_pv_objective(directional::IntrinsicFaceTangentBundle *tb,
                           const Eigen::MatrixXcd &pv_field) {
  int N = pv_field.cols();
  bool sign_symmetry = N % 2 == 0;
  // Number of coefficients per face.
  int nc = sign_symmetry ? N / 2 : N;
  int nf = pv_field.rows();
  int jump = sign_symmetry ? 2 : 1;

  smooth_polyvector_objective->get_quadratic_term(1).update_b(
      pv_field.rows(), [&](int i, const auto &add_b) {
        for (int n = 0; n < nc; n++) {
          add_b(pv_field(i, n));
        }
      });

  smooth_polyvector_objective->get_quadratic_term(2).update_b(
      pv_field.rows(), [&](int i, const auto &add_b) {
        std::complex<double> normalized =
            pv_field(i, 0) / std::abs(pv_field(i, 0));
        add_b(normalized);
      });

  if (const_faces.size() > 1) {
    smooth_polyvector_objective->get_quadratic_term(4).update_b(
        const_faces.size(), [&](int i, const auto &add_b) {
          auto closest = get_closest_aligned_polyvector(tb, pv_field, i);
          for (int n = 0; n < nc; n++) {
            add_b(closest(n));
          }
        });
  }
}

static Eigen::MatrixXcd
smooth_polyvector_field(directional::IntrinsicFaceTangentBundle *tb,
                        const Eigen::MatrixXcd &pv_field) {
  int N = pv_field.cols();
  bool sign_symmetry = N % 2 == 0;
  // Number of coefficients per face.
  int nc = sign_symmetry ? N / 2 : N;
  int nf = pv_field.rows();
  int jump = sign_symmetry ? 2 : 1;

  // Recompute the objective if the parameters have changed.
  if (cached_params != params) {
    compute_smooth_pv_objective(tb, pv_field);
  } else {
    update_smooth_pv_objective(tb, pv_field);
  }

  Eigen::MatrixXcd sol = smooth_polyvector_objective->solve();
  if (!sign_symmetry) {
    return sol;
  }
  Eigen::MatrixXcd res =
      Eigen::MatrixXcd::Zero(pv_field.rows(), pv_field.cols());
  for (int f = 0; f < nf; f++) {
    for (int i = 0; i < N; i += jump) {
      res(f, i) = sol(f, i / jump);
    }
  }

  return res;
}

static void compute_curl_objective(directional::CartesianField &field) {
  directional::IntrinsicFaceTangentBundle *tb =
      (directional::IntrinsicFaceTangentBundle *)field.tb;
  int N = field.N;
  bool has_sign_symmetry = N != 2 && N % 2 == 0;
  int nc = has_sign_symmetry ? N / 2 : N;
  int nf = field.intField.rows();
  int jump = has_sign_symmetry ? 2 : 1;

  reduce_curl_objective = std::make_unique<Optiz::QuadraticObjectiveD>(
      std::vector<int>{nf, nc * 2});
  // Closeness term.
  reduce_curl_objective->add_weighted_equations_with_rhs(
      nf, [&](int i, auto &x, const auto &add_equation) {
        double face_area = tb->mesh->faceAreas(i);
        for (int n = 0; n < nc * 2; n++) {
          add_equation(params.closeness_weight * face_area, x(i, n),
                       field.intField(i, n));
        }
      });

  // Alignment term.
  if (const_faces.size() > 1) {
    reduce_curl_objective->add_weighted_equations_with_rhs(
        const_faces.size(), [&](int i, auto &x, const auto &add_equation) {
          int f = const_faces(i);
          double face_area = tb->mesh->faceAreas(f);
          Eigen::Vector2d vec = get_face_constraint_real(tb, i);
          int closest_match = -1;
          double highest_dot_product = 0;
          for (int n = 0; n < nc; n++) {
            Eigen::Vector2d vf = Eigen::Vector2d(field.intField(f, n * 2),
                                                 field.intField(f, n * 2 + 1))
                                     .normalized();
            double dp = vf.dot(vec);
            if (closest_match == -1 ||
                std::abs(dp) > std::abs(highest_dot_product)) {
              highest_dot_product = dp;
              closest_match = n;
            }
          }
          double sign = highest_dot_product > 0 ? 1 : -1;
          add_equation(params.alignment_weight * face_area,
                       x(f, closest_match * 2), sign * vec(0));
          add_equation(params.alignment_weight * face_area,
                       x(f, closest_match * 2 + 1), sign * vec(1));
        });
  }

  auto start = std::chrono::steady_clock::now();
  // Project curl term.
  reduce_curl_objective->add_hard_constraints(
      tb->adjSpaces.rows(), [&](int i, const auto &x, auto &add_eq) {
        // using T = FACTORY_TYPE(x);
        int f = tb->adjSpaces(i, 0), g = tb->adjSpaces(i, 1);
        if (f == -1 || g == -1)
          return;
        // return Eigen::VectorX<T>(0);

        // Shared edge in local real coordinates.
        auto [ef0, ef1] = get_shared_edge_local_frames_real(tb, i);

        // Term for the curl of each pair of matching vectors between f0, f1.
        // Eigen::VectorX<T> res(nc);
        for (int n = 0; n < nc; n++) {
          // nth vector in face f.
          // Eigen::Vector2<T> vf(x(f, n * 2), x(f, n * 2 + 1));
          Optiz::VariablesVec2 vf(x(f, n * 2), x(f, n * 2 + 1));
          // Matching vector to 'n' in triangle 'g'.
          int m = (n + field.matching(i) + field.N) % field.N;
          // If the matching vector is in the second half of the field, take
          // the one from the first half and flip the sign.
          bool flip = has_sign_symmetry && m >= nc;
          if (flip)
            m -= nc;
          // mth vector in face g.
          // Eigen::Vector2<T> vg(x(g, m * 2), x(g, m * 2 + 1));
          Optiz::VariablesVec2 vg(x(g, m * 2), x(g, m * 2 + 1));
          if (flip) {
            // vg = -vg;
            add_eq(vf.dot(ef0) + vg.dot(ef1));
            // add_eq(x(f, n * 2) * ef0(0) + x(f, n * 2 + 1) * ef0(1) + x(g, m *
            // 2) * ef1(0) + x(g, m * 2 + 1) * ef1(1));
          } else {
            add_eq(vf.dot(ef0) - vg.dot(ef1));
            // add_eq(x(f, n * 2) * ef0(0) + x(f, n * 2 + 1) * ef0(1) - x(g, m *
            // 2) * ef1(0) - x(g, m * 2 + 1) * ef1(1));
          }

          // Enforce zero curl.
          // res(n) = vf.dot(ef0) - vg.dot(ef1) - 0.0;
        }
        // return res;
      });
  auto end = std::chrono::steady_clock::now();
  std::cout << "Curl objective time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
}

static void project_curl(directional::CartesianField &field) {
  compute_curl_objective(field);

  Eigen::MatrixXd res = reduce_curl_objective->solve();
  if (field.N == 2 || field.N % 2 != 0) {
    field.set_intrinsic_field(res);
    return;
  }
  Eigen::MatrixXd int_field(field.intField.rows(), field.intField.cols());
  int_field.leftCols(field.N) = res;
  int_field.rightCols(field.N) = -res;
  field.set_intrinsic_field(int_field);
}

void optimize_rosy_field(directional::CartesianField &field, int max_iters,
                         bool continue_previous) {
  // Init variables.
  stop_optimization = false;
  optimization_updated = false;
  if (!continue_previous) { // Otherwise, continue previous values.
    params = OptimizationParams();
    cached_params.smoothness_weight = -1.0;
  }
  directional::CartesianField tmp_field;
  for (int i = 0; i < max_iters && !stop_optimization; i++) {
    if (params.iter != 0 && params.iter % 5 == 0) {
      params.smoothness_weight *= 0.5;
    }
    // 1. Convert to polyvector field.
    Eigen::MatrixXcd pv_field = raw_to_polyvector(field);
    // 2. Smooth it.
    Eigen::MatrixXcd smoothed_pv_field = smooth_polyvector_field(
        (directional::IntrinsicFaceTangentBundle *)field.tb, pv_field);
    // 3. Convert back to raw field.
    directional::CartesianField directional_pv_field;
    directional_pv_field.init(
        *field.tb, directional::fieldTypeEnum::POLYVECTOR_FIELD, field.N);
    directional_pv_field.set_intrinsic_field(smoothed_pv_field);
    tmp_field = polyvector_to_raw(directional_pv_field, field.N % 2 == 0);

    // 4. Normalize.
    normalize_field(tmp_field);
    // 5. Update matching.
    directional_::principal_matching(tmp_field);

    // 6. Project curl.
    project_curl(tmp_field);
    field = tmp_field;

    optimization_updated = true;
    params.iter++;
    cached_params = params; // Make sure we don't recompute the objective.
  };
  stop_optimization = true;
}

static Eigen::MatrixXcd
real_field_to_complex(const directional::CartesianField &field) {
  Eigen::MatrixXcd res(field.intField.rows(), field.intField.cols() / 2);
  for (int i = 0; i < res.rows(); i++) {
    for (int j = 0; j < res.cols(); j++) {
      res(i, j) = std::complex<double>(field.intField(i, j * 2),
                                       field.intField(i, j * 2 + 1));
    }
  }
  return res;
}

static Eigen::MatrixXd complex_field_to_real(const Eigen::MatrixXcd &field) {
  Eigen::MatrixXd res(field.rows(), field.cols() * 2);
  for (int i = 0; i < field.rows(); i++) {
    for (int j = 0; j < field.cols(); j++) {
      res(i, j * 2) = field(i, j).real();
      res(i, j * 2 + 1) = field(i, j).imag();
    }
  }
  return res;
}

directional::CartesianField
smooth_frame_field(directional::IntrinsicFaceTangentBundle *tb,
                   const directional::CartesianField &field) {
  Eigen::MatrixXcd complex_field = real_field_to_complex(field);
  int nf = complex_field.rows();
  int nc = complex_field.cols();

  smooth_polyvector_objective =
      std::make_unique<Optiz::QuadraticObjectiveCD>(std::vector<int>{nf, nc});
  // Smoothness term.
  smooth_polyvector_objective->add_weighted_equations(
      tb->adjSpaces.rows(), [&](int i, auto &x, const auto &add_equation) {
        int f = tb->adjSpaces(i, 0), g = tb->adjSpaces(i, 1);
        if (f == -1 || g == -1)
          return;
        double weight = params.smoothness_weight * get_edge_weight(tb, i);
        // Get the shared edge in local coords of faces f, g.
        auto [ef, eg] = get_shared_edge_local_frames(tb, i);
        auto efn = std::conj(ef);
        auto egn = std::conj(eg);
        for (int n = 0; n < nc; n++) {
          add_equation(weight, efn * x(f, n) - egn * x(g, n));
        }
      });
  // Closeness term.
  smooth_polyvector_objective->add_weighted_equations_with_rhs(
      nf, [&](int i, auto &x, const auto &add_equation) {
        double weight = params.closeness_weight * tb->mesh->faceAreas(i);
        for (int n = 0; n < nc; n++) {
          add_equation(weight, x(i, n), complex_field(i, n));
        }
      });
  // Orthogonality term.
  auto i = std::complex<double>(0, 1);
  smooth_polyvector_objective->add_weighted_equations(
      nf, [&](int ind, auto &x, const auto &add_equation) {
        double weight = params.rosy_weight * tb->mesh->faceAreas(ind);
        add_equation(weight, x(ind, 1) - i * x(ind, 0));
      });
  Eigen::MatrixXcd sol = smooth_polyvector_objective->solve();

  directional::CartesianField res;
  res.init(*tb, directional::fieldTypeEnum::RAW_FIELD, nc);
  res.set_intrinsic_field(sol);
  return res;
}

void optimize_frame_field(directional::CartesianField &field, int max_iters,
                          bool continue_previous) {
  // Init variables.
  stop_optimization = false;
  optimization_updated = false;
  if (!continue_previous) { // Otherwise, continue previous values.
    params = OptimizationParams();
    cached_params.smoothness_weight = -1.0;
  }
  directional::CartesianField tmp_field;
  for (int i = 0; i < max_iters && !stop_optimization; i++) {
    // if (params.iter != 0 && params.iter % 3 == 0) {
    //   params.smoothness_weight *= 0.5;
    // }
    params.smoothness_weight *= 0.8;
    // 1. Smooth the frame field.
    tmp_field = smooth_frame_field(
        (directional::IntrinsicFaceTangentBundle *)field.tb, field);

    // 2. Normalize.
    normalize_field(tmp_field);

    // 3. Update matching.
    tmp_field.matching = Eigen::VectorXi::Zero(tmp_field.tb->adjSpaces.rows());

    // 4. Project curl.
    project_curl(tmp_field);
    if ((field.intField - tmp_field.intField).norm() / field.intField.norm() <
        1e-3) {
      break;
    }
    field = tmp_field;

    optimization_updated = true;
    params.iter++;
    cached_params = params; // Make sure we don't recompute the objective.
  };
  stop_optimization = true;
}

} // namespace param