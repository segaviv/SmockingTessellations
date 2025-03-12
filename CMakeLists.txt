cmake_minimum_required(VERSION 3.12)
include(FetchContent)

project(smocking_tessellations)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(FETCHCONTENT_QUIET FALSE)

# libigl
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.4.0)
# Directional
FetchContent_Declare(
  Directional
  GIT_REPOSITORY https://github.com/avaxman/Directional.git
  GIT_TAG v2.0.0)
# polyscope
FetchContent_Declare(
  polyscope
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_TAG v2.2.1)
# Optiz
FetchContent_Declare(
  Optiz
  GIT_REPOSITORY https://github.com/segaviv/optiz.git
  GIT_TAG master)
# Eigen3
find_package(Eigen3 NO_MODULE QUIET)
# CGAL.
FetchContent_Declare(
  CGAL URL https://github.com/CGAL/cgal/releases/download/v5.6.2/CGAL-5.6.2.tar.xz
           CONFIGURE_COMMAND "" BUILD_COMMAND "")
FetchContent_Populate(CGAL)
set(CGAL_DIR ${cgal_SOURCE_DIR})
find_package(CGAL REQUIRED)
# ImGuiFileDialog
FetchContent_Declare( 
  ImGuiFileDialog
  GIT_REPOSITORY https://github.com/aiekick/ImGuiFileDialog.git
  GIT_TAG master)

find_package(OpenMP QUIET)

FetchContent_MakeAvailable(libigl polyscope Optiz ImGuiFileDialog Directional)

target_include_directories(ImGuiFileDialog
                           PUBLIC ${polyscope_SOURCE_DIR}/deps/imgui/imgui)

file(
  GLOB
  SRC_FILES
  *.cpp
  parameterization/*.cpp
  smocking_pattern/*.cpp
  libcache/*.cpp
  state/*.cpp
  gui/*.cpp
  view/*.cpp
  utils/*.cpp
  layout/*.cpp
  cut/*.cpp
  io/*.cpp
  visualization/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC ${PROJECT_SOURCE_DIR} ${directional_SOURCE_DIR}/include
         ${directional_SOURCE_DIR}/external/SaddlePoint/include
         ${ImGuiFileDialog_SOURCE_DIR})
if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found, building without it.")
endif()
# ${PROJECT_SOURCE_DIR}
target_link_libraries(
  ${PROJECT_NAME} PUBLIC polyscope igl::core Optiz ImGuiFileDialog 
                         CGAL::CGAL)
