add_executable(sph main.cpp)
target_compile_features(sph
  PRIVATE
    cxx_std_11
  )
target_compile_options(sph
  PRIVATE
    # For GCC and Clang
    $<$<CXX_COMPILER_ID:GNU>:-ffast-math -funroll-loops -Wall -Wextra -Wpedantic>
    # For Visual Studio
    $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /W4>
    # For Intel comipler
    # $<$<CXX_COMPILER_ID:Intel>:$<IF:$<PLATFORM_ID:Windows>,/fp:fast /warn:all,-fp-model fast -warn all>>
  )
target_link_libraries(sph PRIVATE FDPS::FDPS)