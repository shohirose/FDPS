#ifndef FDPS_FUNCTIONS_HPP
#define FDPS_FUNCTIONS_HPP

#include <cstdint>
#include <iostream>

namespace ParticleSimulator {

/// @brief Initialize the library.
///
/// This function calls MPI_Init() and MemoryPool::initialize().
void Initialize(int& argc, char**& argv, int64_t memoryPoolSize = 100'000'000);

/// @brief Finalize the library.
///
/// This function calls MPI_Finalize().
void Finalize();

/// @brief Print license terms
void PrintLicense(std::ostream& os);

}  // namespace ParticleSimulator

#endif  // FDPS_FUNCTIONS_HPP