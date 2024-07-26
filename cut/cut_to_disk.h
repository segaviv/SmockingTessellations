#pragma once

#include <Eigen/Eigen>
#include <set>
#include <vector>

#include <utils/Hmesh.h>
#include <utils/ModifiedMesh.h>
#include <utils/find_duplicates.h>

namespace cut
{
    std::vector<int> singularity_indices(utils::Hmesh& mesh, int symmetry);
    utils::ModifiedMesh to_disk(utils::Hmesh& mesh, int symmetry, bool add_cut_twin_edges = false);
    std::vector<std::set<int>> boundary_loops(utils::Hmesh& mesh, 
    const std::set<int>& additional_boundary_verts);
    std::vector<int> connect_boundaries(utils::Hmesh& mesh,
                                    std::vector<std::set<int>>& boundaries);
} // namespace remesh
