#pragma once

#include "smocking_pattern/UnitSmockingPattern.h"
#include "utils/Hmesh.h"
#include <Eigen/Eigen>
#include <memory>
#include <parameterization/param.h>
#include <smocking_pattern/Tangram.h>
#include <string>
#include <parameterization/pullback.h>
#include <visualization/seamless_arap.h>

namespace state {
// The triangle mesh.
bool load_mesh(const std::string &filename);
extern utils::Hmesh mesh;
extern utils::Hmesh remeshed;
extern utils::Hmesh remeshed_cut;

extern int symmetry;
extern param::FieldData field_data;
extern directional::CartesianField initial_field;

// Pullback result & optimized layout.
extern param::PullbackResult pullback_result;
extern Eigen::MatrixXd optimized_layout;

// ARAP.
extern std::unique_ptr<visualization::SeamlessARAP> pleats_arap;
extern std::unique_ptr<visualization::SeamlessARAP> underlay_arap;
extern Eigen::MatrixXd arap_verts;
extern Eigen::MatrixXi arap_faces;

// Visualizing the planar tangram.
extern int rep_x, rep_y;
extern float stitching_scale;
extern smocking::Tangram tangram;
void set_unit_pattern(smocking::UnitSmockingPattern pattern);

void init_vector_field(int symmetry = 1, bool align_to_boundary = false);

void reset();

} // namespace state