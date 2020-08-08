#pragma once

#include "FDPS/ps_defs.hpp"
#include "FDPS/vector3.hpp"

namespace ParticleSimulator {

class MortonKey {
 private:
  enum {
// MSB is always 0.
// next 3-bits represent octant index of 8 cells with level 1
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    kLevMax = 30,
#else
    kLevMax = 21,
#endif
  };

  F64 half_len_;
  F64vec center_;
  F64 normalized_factor_;

  MortonKey() = default;
  ~MortonKey() = default;
  MortonKey(const MortonKey &) = default;
  MortonKey &operator=(const MortonKey &) = default;

  static MortonKey &getInstance() noexcept;
  static U64 separateBit(const U64 _s_in) noexcept;

 public:
  static void initialize(const F64 half_len,
                         const F64vec &center = 0.0) noexcept;
  static U64 getKey(F64vec pos) noexcept;
  static S32 getCellID(const S32 lev, const U64 mkey) noexcept;
};

}  // namespace ParticleSimulator
