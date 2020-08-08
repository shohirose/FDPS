add_executable(vdwtest vdwtest.cpp)
target_link_libraries(vdwtest PRIVATE FDPS::FDPS)
target_compile_options(vdwtest
  PRIVATE
    # For GCC and Clang
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:-ffast-math -funroll-loops -Wall -Wextra -Wpedantic>
    # For Visual Studio
    $<$<CXX_COMPILER_ID:MSVC>:/fp:fast /W4>
    # For Intel comipler
    # $<$<CXX_COMPILER_ID:Intel>:$<IF:$<PLATFORM_ID:Windows>,/fp:fast /warn:all,-fp-model fast -warn all>>
  )