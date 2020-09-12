#include "FDPS/mpi/communicator.hpp"

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

namespace ParticleSimulator {

namespace mpi {

Communicator::Communicator() : size_{1}, rank_{0} {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
#endif
}

}  // namespace mpi

}  // namespace ParticleSimulator
