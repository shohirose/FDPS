#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
#include <omp.h>
#endif

namespace ParticleSimulator {

class MemoryPool {
 public:
  static size_t getNSegment() noexcept { return getInstance().n_segment_; }
  static size_t getSize() noexcept { return getInstance().size_; }
  static void initialize(const size_t _cap) noexcept;
  static void reInitialize(const size_t size_add) noexcept;
  static void unifyMem() noexcept;
  static void alloc(const size_t _size, int &_id_mpool, void *ptr_data,
                    void *&ret) noexcept;
  static void freeMem(const int id_seg) noexcept;
  static void dump() noexcept;

 private:
  enum {
    ALIGN_SIZE = 8,
    N_SEGMENT_LIMIT = 10000,
  };

  struct EmergencyBuffer {
    void *data;
    size_t cap;
    bool used;
    void *ptr_data;
    void *ptr_id_mpool;
  };

  void *bottom_;
  void *top_;
  size_t cap_;
  size_t size_;
  size_t n_segment_;
  size_t cap_per_seg_[N_SEGMENT_LIMIT];
  bool used_per_seg_[N_SEGMENT_LIMIT];
  void *ptr_data_per_seg_[N_SEGMENT_LIMIT];
  std::vector<EmergencyBuffer> emerg_bufs_;

  MemoryPool() = default;
  ~MemoryPool() = default;
  MemoryPool(const MemoryPool &mem) = default;
  MemoryPool &operator=(const MemoryPool &mem) = default;

  static MemoryPool &getInstance() noexcept;
  static size_t getAlignSize(const size_t _size) noexcept;
  static bool isLastSegment(const int id_seg) noexcept;
  static bool inParallelRegion() noexcept;
};

}  // namespace ParticleSimulator
