#ifndef FDPS_MPI_COMMUNICATOR_HPP
#define FDPS_MPI_COMMUNICATOR_HPP

namespace ParticleSimulator {

namespace mpi {

/// @brief MPI communicator
///
/// mpi::Environment must be instantiated in advance.
class Communicator {
 public:
  Communicator();

  /// @brief Number of MPI processes
  int size() const noexcept { return size_; }

  /// @brief The rank of the current MPI process
  int rank() const noexcept { return rank_; }

 private:
  int size_;
  int rank_;
};

}  // namespace mpi

}  // namespace ParticleSimulator

#endif  // FDPS_MPI_COMMUNICATOR_HPP