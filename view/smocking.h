#pragma once
#include <smocking_pattern/UnitSmockingPattern.h>
#include <smocking_pattern/Tangram.h>

namespace view {

void visualize_smocking_pattern(const smocking::UnitSmockingPattern& pattern);

void visualize_tangram(const smocking::Tangram& tangram, int repeat_x = 1, int repeat_y = 1);

void show_pulled_back_tangram();
void show_open_configuration();
void show_optimized_layout();
void show_stitching_lines();
void show_arap_result();

}