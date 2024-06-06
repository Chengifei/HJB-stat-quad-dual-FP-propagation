cmake_minimum_required(VERSION 3.10)
include(FindPkgConfig)
project(HJB_solve)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

pkg_check_modules(PKG REQUIRED eigen3)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

include_directories(SYSTEM ${PKG_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
link_directories(${PKG_LIBRARY_DIRS})
link_libraries(${PKG_LIBRARIES})
add_compile_definitions(NPY_NO_DEPRECATED_API=NPY_1_20_API_VERSION EIGEN_DONT_PARALLELIZE EIGEN_MAX_ALIGN_BYTES=0 NDEBUG)

add_library(HJB_solve SHARED N_tilde.cpp wrapper.cpp)
set_target_properties(HJB_solve PROPERTIES SUFFIX .so PREFIX "")
target_compile_options(HJB_solve PRIVATE -march=native)
