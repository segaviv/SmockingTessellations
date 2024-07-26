#include "main_menu.h"
#include "cut/auxetic_pattern.h"
#include "cut/cut_boundary_faces.h"
#include "cut/cut_to_disk.h"
#include "imgui.h"
#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/screenshot.h"
#include "polyscope/view.h"
#include "smocking_pattern/UnitSmockingPattern.h"
#include "utils/Hmesh.h"
#include "utils/conversions.h"
#include "utils/slicing.h"
#include <ImGuiFileDialog.h>
#include <igl/upsample.h>
#include <io/save_mesh.h>
#include <layout/optimize.h>
#include <memory>
#include <parameterization/grid_param.h>
#include <parameterization/param.h>
#include <parameterization/pullback.h>
#include <parameterization/remeshing.h>
#include <smocking_pattern/Tangram.h>
#include <state/state.h>
#include <thread>
#include <utils/find_duplicates.h>
#include <view/manager.h>
#include <view/smocking.h>
#include <visualization/prepare_arap.h>
#include <visualization/seamless_arap.h>

namespace gui {

static const char *patterns[] = {"Arrow", "Leaf", "braid", "box",
                                 "twisted_square", "heart", "brick"};
static int selected_pattern = 0;
static float rotate_uv_angle = 0.0;
static float scale_uv = 1.0;
static Eigen::Vector2f shift_uv = Eigen::Vector2f::Zero();
static const char *symmetries[] = {"None", "3", "4", "6"};
static const int symmetry_values[] = {1, 3, 4, 6};
static int selected_symmetry = 0;
static bool align_to_boundary = false;

#define DISABLED_BUTTON(name, cond, code)                                      \
  ImGui::BeginDisabled(cond);                                                  \
  if (ImGui::Button(name)) {                                                   \
    code                                                                       \
  }                                                                            \
  ImGui::EndDisabled();

static void open_mesh_dialog() {
  if (ImGui::Button("Load mesh")) {
    IGFD::FileDialogConfig config{.flags = ImGuiFileDialogFlags_Modal};
    ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose mesh",
                                            ".obj,.off", config);
  }
  if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey")) {
    if (ImGuiFileDialog::Instance()->IsOk()) {
      std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
      state::load_mesh(filePathName);
      polyscope::removeAllStructures();
      view::add_surface_mesh(state::mesh);
    }
    ImGuiFileDialog::Instance()->Close();
  }
}

static void parameterization_menu() {
  if (ImGui::CollapsingHeader("[Debug] parameterization")) {

    if (ImGui::BeginCombo("Rotational symmetry", symmetries[selected_symmetry],
                          ImGuiComboFlags_WidthFitPreview)) {
      for (int i = 0; i < IM_ARRAYSIZE(symmetries); i++) {
        if (ImGui::Selectable(symmetries[i], i == selected_symmetry)) {
          selected_symmetry = i;
          state::symmetry = symmetry_values[selected_symmetry];
        }
      }
      ImGui::EndCombo();
    }
    // Init vector field.
    if (ImGui::Button("Init")) {
      state::init_vector_field(symmetry_values[selected_symmetry],
                               align_to_boundary);
      view::show_vector_field(state::field_data);
    }
    ImGui::SameLine();
    ImGui::Checkbox("Align to boundary", &align_to_boundary);

    // Optimize.
    DISABLED_BUTTON("Optimize field", !param::stop_optimization, {
      std::thread([&] {
        if (state::field_data.field.N > 2) {
          param::optimize_rosy_field(state::field_data.field,
                                     100, true);
        } else {
          param::optimize_frame_field(state::field_data.field, 100, true);
        }
      }).detach();
    });
    ImGui::SameLine();
    DISABLED_BUTTON("Stop", param::stop_optimization,
                    { param::stop_optimization = true; });

    // Integrate and remesh.
    if (ImGui::Button("Remesh")) {
      param::setup_integration(state::field_data.field);
      param::integrate(state::field_data.field);
      param::mesh_isolines(state::field_data.mesh);
      state::remeshed =
          utils::Hmesh(param::VPolyMesh, param::FPolyMesh, param::DPolyMesh);
      if (state::field_data.field.N == 6) {
        view::add_surface_mesh(state::remeshed.dual(), "remeshed");
        view::show_only("remeshed");
      } else {
        view::add_surface_mesh(state::remeshed, "remeshed");
        view::show_only("remeshed");
      }
    }
    if (ImGui::Button("Show singularities")) {
      auto inds = slice::indices(state::remeshed.nv(), [&](int i) {
        return !state::remeshed.is_boundary_vertex(i) &&
               state::remeshed.verts_edges[i].size() != 6;
      });
      polyscope::registerPointCloud("singularities",
                                    state::remeshed.V(inds, Eigen::all));
    }
    if (ImGui::Button("cut to disk")) {
      auto cut_mesh = cut::to_disk(state::remeshed, state::symmetry);
      state::remeshed_cut = cut_mesh.mesh;
      auto dups = utils::get_duplicated_verts(state::remeshed_cut);
      std::vector<int> inds;
      for (auto [a, b] : dups) {
        inds.push_back(a);
      }
      polyscope::registerPointCloud("cut",
                                    state::remeshed_cut.V(inds, Eigen::all));
      state::remeshed_cut = cut::boundary_faces(state::remeshed_cut);
      state::pullback_result = cut::cut_auxetic_pattern(state::remeshed_cut);
      view::show_pulled_back_tangram();
    }
    ImGui::SameLine();
    ImGui::InputDouble("length ratio", &param::length_ratio);
  }
}

static smocking::UnitSmockingPattern &get_selected_pattern() {
  switch (selected_pattern) {
  case 0:
    return smocking::UnitSmockingPattern::ARROW;
  case 1:
    return smocking::UnitSmockingPattern::LEAF;
  case 2:
    return smocking::UnitSmockingPattern::BRAID;
  case 3:
    return smocking::UnitSmockingPattern::BOX;
  case 4:
    return smocking::UnitSmockingPattern::TWISTED_SQUARE;
  case 5:
    return smocking::UnitSmockingPattern::HEART;
  case 6:
    return smocking::UnitSmockingPattern::BRICK;
  default:
    return smocking::UnitSmockingPattern::ARROW;
  }
}

static void export_mesh(const Eigen::MatrixXd &verts,
                        const Eigen::MatrixXi &faces, const Eigen::MatrixXd &uv,
                        const std::string &filename) {
  std::fstream s{filename + ".obj", s.binary | s.trunc | s.in | s.out};

  for (int i = 0; i < verts.rows(); i++) {
    s << "v " << verts(i, 0) << " " << verts(i, 1) << " " << verts(i, 2)
      << std::endl;
  }
  for (int i = 0; i < uv.rows(); i++) {
    s << "vt " << uv(i, 0) << " " << uv(i, 1) << std::endl;
  }
  for (int i = 0; i < faces.rows(); i++) {
    s << "f " << faces(i, 0) + 1 << "/" << faces(i, 0) + 1 << " "
      << faces(i, 1) + 1 << "/" << faces(i, 1) + 1 << " " << faces(i, 2) + 1
      << "/" << faces(i, 2) + 1 << " " << std::endl;
  }
  s.close();
}

static void visualization_menu() {
  ImGui::PushItemWidth(100);
  if (ImGui::CollapsingHeader("Visualization",
                              ImGuiTreeNodeFlags_DefaultOpen)) {
    static bool show_vector_field = false;
    auto show_original = [] {
      view::remove_all();
      view::add_surface_mesh(state::mesh, "mesh");
    };
    auto show_vf = [] {
      if (state::field_data.field.intField.size())
        view::show_vector_field(state::field_data, "mesh");
    };
    if (ImGui::Button("Show original mesh")) {
      show_original();
      if (show_vector_field)
        show_vf();
    }
    ImGui::SameLine();
    ImGui::Checkbox("Show vector field", &show_vector_field);
    DISABLED_BUTTON("Show closed configuration",
                    state::pullback_result.vertices.empty(),
                    { view::show_pulled_back_tangram(); });
    ImGui::SameLine();
    DISABLED_BUTTON("Show open configuration",
                    state::pullback_result.vertices.empty(),
                    { view::show_open_configuration(); });
    DISABLED_BUTTON("Show optimized layout",
                    state::optimized_layout.size() == 0,
                    { view::show_optimized_layout(); });
    ImGui::SameLine();
    DISABLED_BUTTON("Show stitching lines", state::optimized_layout.size() == 0,
                    { view::show_stitching_lines(); });
    DISABLED_BUTTON("Export ARAP result", state::arap_verts.size() == 0, {
      export_mesh(state::arap_verts,
                  convert::to_eig_mat(state::pleats_arap->faces()),
                  state::arap_verts, "arap_result.obj");
    });
  }
  ImGui::PopItemWidth();
}

static void smocking_pattern_menu() {
  ImGui::PushItemWidth(100);
  if (ImGui::CollapsingHeader("Debug Smocking pattern")) {

    if (ImGui::BeginCombo("Smocking pattern", patterns[selected_pattern],
                          ImGuiComboFlags_WidthFitPreview)) {
      bool changed = false;
      for (int i = 0; i < IM_ARRAYSIZE(patterns); i++) {
        if (ImGui::Selectable(patterns[i], i == selected_pattern)) {
          if (selected_pattern != i) {
            changed = true;
            selected_pattern = i;
          }
          break;
        }
      }
      if (changed) {
        state::set_unit_pattern(get_selected_pattern());
        state::tangram.set_stitching_scale(state::stitching_scale);
        view::visualize_tangram(state::tangram, state::rep_x, state::rep_y);
      }
      ImGui::EndCombo();
    }
    ImGui::SameLine();
    if (ImGui::SliderFloat("Stitching Scale", &state::stitching_scale, 0.,
                           1.0)) {
      state::stitching_scale = std::max(state::stitching_scale, 1e-4f);
      state::tangram.set_stitching_scale(state::stitching_scale);
      view::visualize_tangram(state::tangram, state::rep_x, state::rep_y);
    }
    if (ImGui::SliderInt("Repeat X", &state::rep_x, 1, 10)) {
      view::visualize_tangram(state::tangram, state::rep_x, state::rep_y);
    }
    ImGui::SameLine();
    if (ImGui::SliderInt("Repeat Y", &state::rep_y, 1, 10)) {
      view::visualize_tangram(state::tangram, state::rep_x, state::rep_y);
    }
  }
  ImGui::PopItemWidth();
}

static void quick_smocking() {
  if (ImGui::BeginCombo("Symmetry", symmetries[selected_symmetry],
                        ImGuiComboFlags_WidthFitPreview)) {
    for (int i = 0; i < IM_ARRAYSIZE(symmetries); i++) {
      if (ImGui::Selectable(symmetries[i], i == selected_symmetry)) {
        selected_symmetry = i;
        state::symmetry = symmetry_values[selected_symmetry];
      }
    }
    ImGui::EndCombo();
  }
  ImGui::SameLine();
  // Add selection of non-symmetric patterns.
  if (selected_symmetry == 0) {
    if (ImGui::BeginCombo("Pattern", patterns[selected_pattern],
                          ImGuiComboFlags_WidthFitPreview)) {
      for (int i = 0; i < IM_ARRAYSIZE(patterns); i++) {
        if (ImGui::Selectable(patterns[i], i == selected_pattern)) {
          selected_pattern = i;
          break;
        }
      }
      state::set_unit_pattern(get_selected_pattern());
      state::tangram.set_stitching_scale(state::stitching_scale);
      ImGui::EndCombo();
    }
  } else {
    ImGui::Checkbox("Align boundary", &align_to_boundary);
  }

  auto pull_back = [&] {
    if (selected_symmetry == 0) {
      if (state::field_data.field.intField.size() == 0) {
        state::init_vector_field(1, false);
        param::optimize_frame_field(state::field_data.field, 50, true);
        param::optimization_updated = false;
      }
      if (state::field_data.NFunction.size() == 0) {
        param::setup_integration(state::field_data.field);
        param::integrate(state::field_data.field);
      }
      param::pullback_no_rosy(rotate_uv_angle, scale_uv,
                              shift_uv.cast<double>().eval());
      view::show_pulled_back_tangram();
    } else if (selected_symmetry == 1) {
      // Optimize the vector field.
      state::init_vector_field(symmetry_values[selected_symmetry],
                               align_to_boundary);
      param::optimize_rosy_field(state::field_data.field, 100, false);
      param::optimization_updated = false;
      // Integrate the vector field.
      param::setup_integration(state::field_data.field);
      param::integrate(state::field_data.field);
      // Remesh.
      param::mesh_isolines(state::field_data.mesh);
      state::remeshed =
          utils::Hmesh(param::VPolyMesh, param::FPolyMesh, param::DPolyMesh);

      // Cut pattern.
      auto cut_mesh = cut::to_disk(state::remeshed, state::symmetry);
      state::remeshed_cut = cut_mesh.mesh;
      state::remeshed_cut = cut::boundary_faces(state::remeshed_cut);
      state::pullback_result = cut::cut_auxetic_pattern(state::remeshed_cut);
      view::show_pulled_back_tangram();
    } else {
      polyscope::warning("Only implemented for 3 symmetries.");
    }
  };
  if (ImGui::Button("Pullback pattern")) {
    pull_back();
  }
  if (selected_symmetry != 0) {
    ImGui::SameLine();
    ImGui::InputDouble("Scale", &param::length_ratio);
  }

  if (selected_symmetry == 0) {
    ImGui::PushItemWidth(100);
    if (ImGui::SliderFloat("Rotation", &rotate_uv_angle, 0, 2 * M_PI)) {
      if (state::field_data.field.intField.size())
        pull_back();
    }
    ImGui::SameLine();
    if (ImGui::SliderFloat("scale", &scale_uv, 0, 10)) {
      if (state::field_data.field.intField.size())
        pull_back();
    }
    ImGui::SameLine();
    if (ImGui::SliderFloat2("shift", shift_uv.data(), -1, 1)) {
      if (state::field_data.field.intField.size())
        pull_back();
    }
    ImGui::PopItemWidth();
  }
  DISABLED_BUTTON(
      "Optimize layout",
      !layout::stop_optimization || !state::pullback_result.vertices.size(),
      { std::thread([&] { layout::optimize_for_surface(); }).detach(); });
  ImGui::SameLine();
  DISABLED_BUTTON("Stop layout optimization", layout::stop_optimization,
                  { layout::stop_optimization = true; });
  static int iter = 0;
  DISABLED_BUTTON("Prepare ARAP",
                  !state::optimized_layout.size() || !layout::stop_optimization,
                  {
                    visualization::prepare_arap();
                    iter = 0;
                  });
  DISABLED_BUTTON("Run ARAP", !state::pleats_arap, {
    if (iter++ == 0) {
      state::arap_verts =
          state::pleats_arap->solve(state::pleats_arap->verts(), true);
    } else {
      state::arap_verts = state::underlay_arap->solve(state::arap_verts);
    }
    view::show_arap_result();
  });
}

bool show_main_menu = true;
void main_menu() {
  open_mesh_dialog();
  quick_smocking();
  ImGui::Separator();
  visualization_menu();
  ImGui::Separator();
  parameterization_menu();
  smocking_pattern_menu();
}
} // namespace gui