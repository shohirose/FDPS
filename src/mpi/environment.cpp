#include "FDPS/mpi/environment.hpp"

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

namespace ParticleSimulator {

namespace mpi {

Environment::Environment([[maybe_unused]] int argc,
                         [[maybe_unused]] char** argv) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Init(&argc, &argv);
#endif
}

Environment::~Environment() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Finalize();
#endif
}

}  // namespace mpi

}  // namespace ParticleSimulator
