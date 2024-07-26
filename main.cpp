#include "imgui.h"
#include "polyscope/options.h"
#include "polyscope/surface_mesh.h"
#include "utils/Hmesh.h"
#include <parameterization/param.h>
#include <polyscope/polyscope.h>
#include <smocking_pattern/Tangram.h>
#include <smocking_pattern/UnitSmockingPattern.h>

#include <gui/main_menu.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <state/state.h>
#include <utils/slicing.h>
#include <view/manager.h>
#include <view/smocking.h>
// #include <Optiz/Common/SparseUtils.h>

void callback() {
  gui::main_menu();
  view::update_visualization();
}

void init_polyscope() {
  polyscope::init();
  polyscope::options::groundPlaneEnabled = false;
  polyscope::options::ssaaFactor = 4;
  polyscope::state::userCallback = callback;
}

int main(int argc, char **argv) {
  init_polyscope();
  if (argc > 1) {
    std::string filename(argv[1]);
    state::load_mesh(filename);
    view::add_surface_mesh(state::mesh);
  }
  Eigen::SparseMatrix<double> bla;

  std::vector<Eigen::Vector2d> test(5);
  for (int i = 0; i <5; i++) {
    test[i] = Eigen::Vector2d(i, i+1);
  }
  std::cout << Eigen::Matrix<double, -1, 2, Eigen::RowMajor>::Map((double*)test.data(), test.size(), 2) << std::endl;

  polyscope::show();
  return 0;
}