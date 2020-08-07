cmake_minimum_required(VERSION 3.8.2)
project(FDPS CXX C)

add_library(FDPS INTERFACE)
target_include_directories(FDPS INTERFACE src)
target_compile_features(FDPS
  INTERFACE
    # C++11
    # cxx_std_11
    # C++14
    # cxx_std_14
    # C++17
    cxx_std_17
  )
add_subdirectory(src)

add_library(FDPS::FDPS ALIAS FDPS)

add_subdirectory(sample/c++/sph)