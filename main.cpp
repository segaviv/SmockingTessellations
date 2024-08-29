#include "polyscope/options.h"
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
  polyscope::show();
  return 0;
}