#ifndef FDPS_MPI_ENVIRONMENT_HPP
#define FDPS_MPI_ENVIRONMENT_HPP

namespace ParticleSimulator {

namespace mpi {

/// @brief Initialize and finalize the MPI environment.
class Environment {
 public:
  Environment(int argc, char** argv);
  ~Environment();
};

}  // namespace mpi

}  // namespace ParticleSimulator

#endif  // FDPS_MPI_ENVIRONMENT_HPP
