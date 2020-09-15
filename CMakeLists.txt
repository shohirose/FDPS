cmake_minimum_required(VERSION 3.8.2)
project(FDPS CXX)

cmake_policy(SET CMP0076 NEW)

find_package(MPI)

add_library(FDPS)
target_include_directories(FDPS PUBLIC include)
target_compile_features(FDPS
  PUBLIC
    # C++11
    # cxx_std_11
    # C++14
    # cxx_std_14
    # C++17
    cxx_std_17
  )
target_compile_definitions(FDPS
  PUBLIC
    $<$<BOOL:${MPI_FOUND}>:PARTICLE_SIMULATOR_MPI_PARALLEL>
  )
target_link_libraries(FDPS
  PUBLIC
    $<$<BOOL:${MPI_FOUND}>:MPI::MPI_CXX>
  )
add_subdirectory(include/FDPS)
add_subdirectory(src)

add_library(FDPS::FDPS ALIAS FDPS)

add_subdirectory(sample/c++/sph)
# add_subdirectory(sample/c++/vdw-test)
# add_subdirectory(src/particle_mesh)