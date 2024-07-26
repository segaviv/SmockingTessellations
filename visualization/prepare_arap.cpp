#include "prepare_arap.h"
#include "utils/conversions.h"

#include <igl/upsample.h>
#include <state/state.h>
#include <unordered_map>
#include <utils/construct.h>
#include <utils/Hmesh.h>

namespace visualization {

std::unordered_map<int, int>
find_duplicated_verts(const utils::Hmesh &mesh, const utils::Hmesh &arap_mesh) {
  std::set<int> covered_verts;
  auto dup = [&](int v) {
    auto it = state::pullback_result.duplicate_verts.find(v);
    if (it == state::pullback_result.duplicate_verts.end())
      return -1;
    return *it->second.begin();
  };
  auto next_bv = [&](int v) { return mesh.verts[v].edge()->next()->vi; };
  auto arap_next_bv = [&](int v) {
    return arap_mesh.verts[v].edge()->next()->vi;
  };
  auto walk_boundary = [&](int start, int end) {
    std::vector<int> path;
    while (start != end) {
      path.push_back(start);
      start = arap_next_bv(start);
    }
    path.push_back(end);
    return path;
  };

  std::unordered_map<int, int> result;
  for (auto &[v, d] : state::pullback_result.duplicate_verts) {
    if (!covered_verts.insert(v).second)
      continue;
    int dup_v = *d.begin();
    result[v] = dup_v;

    int next_v = next_bv(v);
    int dup_next_v = dup(next_v);
    if (dup_next_v == -1) {
      if (next_bv(next_v) == dup_v)
        dup_next_v = next_v;
      else {
        continue;
      }
    }
    if (next_bv(dup_next_v) != dup_v) {
      std::cout << "ERROR: next_bv(dup_next_v) != dup_v" << std::endl;
      continue;
    }

    // v -> next_v... dup(v) -> dup_next_v.
    std::vector<int> verts = walk_boundary(v, next_v);
    std::vector<int> dup_verts = walk_boundary(dup_next_v, dup_v);
    if (verts.size() != dup_verts.size()) {
      std::cout << "ERROR: verts.size() != dup_verts.size()" << std::endl;
      continue;
    }

    for (int i = 0; i < verts.size(); i++) {
      result[verts[i]] = dup_verts[verts.size() - i - 1];
      result[dup_verts[verts.size() - i - 1]] = verts[i];
    }
    covered_verts.insert(dup_next_v);
  }

  return result;
}

void prepare_arap() {
  Eigen::MatrixXd V = convert::to_3d(
      state::optimized_layout(state::pullback_result.inside_verts_indices,
                              Eigen::all)
          .eval());
  Eigen::MatrixXi F = convert::to_eig_mat(state::pullback_result.reduced_faces);
  utils::Hmesh mesh(V, F);
  igl::upsample(V, F, 3);
  // ARAP mesh.
  utils::Hmesh arap_mesh(V, F);
  // Update duplicate vertices.
  auto dup = find_duplicated_verts(mesh, arap_mesh);

  Eigen::MatrixXd closed_config =
      state::pullback_result.get_closed_configuration();
  closed_config =
      closed_config(state::pullback_result.inside_verts_indices, Eigen::all);

  // Create arap.
  state::pleats_arap = std::make_unique<visualization::SeamlessARAP>(arap_mesh, dup);
  state::pleats_arap->set_objective();
  // Same objective as pleat ARAP.

  state::underlay_arap =
      std::make_unique<visualization::SeamlessARAP>(*state::pleats_arap.get());

  state::pleats_arap->set_fixed(construct::linspace<int>(closed_config.rows()),
                                closed_config);
  state::pleats_arap->precompute();

  state::underlay_arap->set_fixed(
      state::pullback_result.inside_underlay_indices,
      closed_config(state::pullback_result.inside_underlay_indices,
                    Eigen::all));
  state::underlay_arap->precompute();
}
} // namespace visualization