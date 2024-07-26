#pragma once

#include "directional/IntrinsicFaceTangentBundle.h"
#include "directional/TriMesh.h"
#include <directional/CartesianField.h>

namespace param {


struct FieldData {
  directional::CartesianField field;
  directional::IntrinsicFaceTangentBundle tangent_bundle;
  directional::TriMesh mesh;
  Eigen::MatrixXd NFunction, NCornerFunction;
};

extern bool stop_optimization;
extern bool optimization_updated;

/**
 * @brief Initializes the field.
 *
 * @param V vertices of the mesh.
 * @param F faces of the mesh.
 * @param symmetry pattern symmetry (1 - no symmetry, otherwise should be
 * 3/4/6).
 * @param const_faces indices of faces with fixed vectors.
 * @param const_vectors fixed vectors for the faces in const_faces.
 */
void init_field(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                FieldData &res, int symmetry = 1,
                bool align_to_boundary = false);

void optimize_rosy_field(directional::CartesianField &field, int max_iters = 100,
                    bool continue_previous = false);

void optimize_frame_field(directional::CartesianField &field, int max_iters = 100,
                    bool continue_previous = false);

} // namespace param