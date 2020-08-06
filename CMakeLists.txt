cmake_minimum_required(VERSION 3.8.2)
project(FDPS CXX C)

add_library(FDPS INTERFACE)
target_include_directories(FDPS INTERFACE src)
target_compile_features(FDPS
  INTERFACE
    cxx_std_11
  )

add_library(FDPS::FDPS ALIAS FDPS)

add_subdirectory(sample/c++/sph)