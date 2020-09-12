#include "FDPS/comm.hpp"

namespace ParticleSimulator {

Comm::Comm() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc_);
#else
  rank_ = 0;
  n_proc_ = 1;
#endif
  //#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
  n_thread_ = omp_get_max_threads();
#else
  n_thread_ = 1;
#endif
}

Comm & Comm::getInstance() {
  static Comm inst;
  return inst;
}

}  // namespace ParticleSimulator
