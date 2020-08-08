add_executable(sph main.cpp)
target_compile_options(sph
  PRIVATE
    # For GCC and Clang
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-ffast-math -funroll-loops -Wall -Wextra -Wpedantic>
    # For Visual Studio
    $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /W4>
    # For Intel comipler
    # $<$<CXX_COMPILER_ID:Intel>:$<IF:$<PLATFORM_ID:Windows>,/fp:fast /warn:all,-fp-model fast -warn all>>
  )
target_link_libraries(sph PRIVATE FDPS::FDPS)