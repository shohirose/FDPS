#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <functional>
#include <map>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include <time.h>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#endif

#include "FDPS/matrix2.hpp"
#include "FDPS/matrix_sym2.hpp"
#include "FDPS/matrix_sym3.hpp"
#include "FDPS/memory_pool.hpp"
#include "FDPS/orthotope2.hpp"
#include "FDPS/orthotope2i.hpp"
#include "FDPS/orthotope3.hpp"
#include "FDPS/orthotope3i.hpp"
#include "FDPS/vector2.hpp"
#include "FDPS/vector3.hpp"

#define PS_DEBUG_CALL(func)                                                   \
  do {                                                                        \
    try {                                                                     \
      func;                                                                   \
      std::cout << "[FDPS msg] " << #func << ": " << getMemSizeUsed() << ", " \
                << Comm::getRank() << ", " << typeid(TSM).name() << ". "      \
                << std::endl;                                                 \
    } catch (std::bad_alloc & e) {                                            \
      std::cout << "[FDPS error] " << #func << ": " << getMemSizeUsed()       \
                << ", " << Comm::getRank() << ", " << typeid(TSM).name()      \
                << ". " << std::endl;                                         \
      MPI_Abort(MPI_COMM_WORLD, 9);                                           \
      std::exit(1);                                                           \
    } catch (...) {                                                           \
      std::cout << "[FDPS unknown error] " << #func << ": "                   \
                << getMemSizeUsed() << ", " << Comm::getRank() << ", "        \
                << typeid(TSM).name() << ". " << std::endl;                   \
      MPI_Abort(MPI_COMM_WORLD, 9);                                           \
      std::exit(1);                                                           \
    }                                                                         \
    MPI_Barrier(MPI_COMM_WORLD);                                              \
    if (Comm::getRank() == 0) std::cout << #func << " passed." << std::endl;  \
  } while (0);

#define PARTICLE_SIMULATOR_PRINT_ERROR(msg)                             \
  {                                                                     \
    std::cout << "PS_ERROR: " << msg << " \n"                           \
              << "function: " << __FUNCTION__ << ", line: " << __LINE__ \
              << ", file: " << __FILE__ << std::endl;                   \
  }

#define PARTICLE_SIMULATOR_PRINT_LINE_INFO()                            \
  {                                                                     \
    std::cout << "function: " << __FUNCTION__ << ", line: " << __LINE__ \
              << ", file: " << __FILE__ << std::endl;                   \
  }

namespace ParticleSimulator {
static const long long int LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll << 30;
static inline void Abort(const int err = -1) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Abort(MPI_COMM_WORLD, err);
#else
  exit(err);
#endif
}
}  // namespace ParticleSimulator

#include "FDPS/reallocatable_array.hpp"

namespace ParticleSimulator {
typedef int S32;
typedef unsigned int U32;
#ifdef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
typedef double F32;
// typedef long long int S32;
// typedef unsigned long long int U32;
#else
typedef float F32;
// typedef int S32;
// typedef unsigned int U32;
#endif
typedef long long int S64;
typedef unsigned long long int U64;
typedef double F64;
typedef Vector2<S32> S32vec2;
typedef Vector3<S32> S32vec3;
typedef Vector2<U64> U64vec2;
typedef Vector3<U64> U64vec3;
typedef Vector2<F32> F32vec2;
typedef Vector3<F32> F32vec3;
typedef Vector2<F64> F64vec2;
typedef Vector3<F64> F64vec3;
typedef MatrixSym2<F32> F32mat2;
typedef MatrixSym3<F32> F32mat3;
typedef MatrixSym2<F64> F64mat2;
typedef MatrixSym3<F64> F64mat3;
typedef Orthotope2i<S32> S32ort2;
typedef Orthotope3i<S32> S32ort3;
typedef Orthotope2<F32> F32ort2;
typedef Orthotope3<F32> F32ort3;
typedef Orthotope2<F64> F64ort2;
typedef Orthotope3<F64> F64ort3;

static const S32 DIMENSION_LIMIT = 3;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
typedef S32vec2 S32vec;
typedef F32vec2 F32vec;
typedef F64vec2 F64vec;
typedef F32mat2 F32mat;
typedef F64mat2 F64mat;
typedef S32ort2 S32ort;
typedef F32ort2 F32ort;
typedef F64ort2 F64ort;
static const S32 DIMENSION = 2;
static const S32 N_CHILDREN = 4;
static const S32 TREE_LEVEL_LIMIT = 30;
static const F32vec SHIFT_CENTER[N_CHILDREN] = {
    F32vec(-0.5, -0.5), F32vec(-0.5, 0.5), F32vec(0.5, -0.5), F32vec(0.5, 0.5)};
#else
typedef S32vec3 S32vec;
typedef F32vec3 F32vec;
typedef F64vec3 F64vec;
typedef F32mat3 F32mat;
typedef F64mat3 F64mat;
typedef S32ort3 S32ort;
typedef F32ort3 F32ort;
typedef F64ort3 F64ort;
static const S32 DIMENSION = 3;
static const S32 N_CHILDREN = 8;
static const S32 TREE_LEVEL_LIMIT = 21;
static const F32vec SHIFT_CENTER[N_CHILDREN] = {
    F32vec(-0.5, -0.5, -0.5), F32vec(-0.5, -0.5, 0.5), F32vec(-0.5, 0.5, -0.5),
    F32vec(-0.5, 0.5, 0.5),   F32vec(0.5, -0.5, -0.5), F32vec(0.5, -0.5, 0.5),
    F32vec(0.5, 0.5, -0.5),   F32vec(0.5, 0.5, 0.5)};
#endif
#ifdef PARTICLE_SIMULATOR_SPMOM_F32
typedef S32 SSP;
typedef F32 FSP;
typedef F32vec FSPvec;
typedef F32mat FSPmat;
#else
typedef S64 SSP;
typedef F64 FSP;
typedef F64vec FSPvec;
typedef F64mat FSPmat;
#endif
typedef U64 CountT;

static const F64 LARGE_DOUBLE = std::numeric_limits<F64>::max() * 0.0625;
static const F64 LARGE_FLOAT = std::numeric_limits<F32>::max() * 0.0625;
static const S64 LARGE_INT = std::numeric_limits<S32>::max() * 0.0625;

///// A.Tanikawa modified from
// In the upper line, the right-hand side is interpreted to be a 32bit-integer.
// static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1<<31;
// static const S64 LIMIT_NUMBER_OF_TREE_PARTICLE_PER_NODE = 1ll<<30;
///// A.Tanikawa modified to

//////////////////
/// enum
enum INTERACTION_LIST_MODE {
  MAKE_LIST,
  MAKE_LIST_FOR_REUSE,
  REUSE_LIST,
};

enum SEARCH_MODE {
  LONG_NO_CUTOFF,
  LONG_CUTOFF,
  LONG_SCATTER,         // new for P^3T
  LONG_CUTOFF_SCATTER,  // new for P^3T + PM
  LONG_SYMMETRY,
  SHORT_GATHER,
  SHORT_SCATTER,
  SHORT_SYMMETRY,
};

enum FORCE_TYPE {
  FORCE_TYPE_LONG,
  FORCE_TYPE_SHORT,
};
struct TagForceLong {
  enum {
    force_type = FORCE_TYPE_LONG,
  };
};
struct TagForceShort {
  enum {
    force_type = FORCE_TYPE_SHORT,
  };
};

struct TagSearchBoundaryConditionOpenOnly {};
struct TagSearchBoundaryConditionOpenPeriodic {};

struct TagSearchLong {};
struct TagSearchLongCutoff {};
struct TagSearchLongScatter {};
struct TagSearchLongSymmetry {};
struct TagSearchLongCutoffScatter {};
struct TagSearchShortGather {};
struct TagSearchShortScatter {};
struct TagSearchShortSymmetry {};

struct TagIpgLongNormal {};
struct TagIpgIn {};
struct TagIpgInAndOut {};
struct TagIpgOut {};

struct TagWithoutCutoff {};
struct TagWithCutoff {};

struct TagNeighborSearchSymmetry {};
struct TagNeighborSearchGather {};
struct TagNeighborSearchScatter {};
struct TagNeighborSearchNo {};

struct SEARCH_MODE_LONG {
  typedef TagForceLong force_type;
  typedef TagSearchLong search_type;
  typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
  typedef TagIpgLongNormal ipg_type;
  typedef TagNeighborSearchNo neighbor_search_type;
  enum {
    search_type_id = LONG_NO_CUTOFF,
  };
};
struct SEARCH_MODE_LONG_CUTOFF {
  typedef TagForceLong force_type;
  typedef TagSearchLongCutoff search_type;
  typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
  typedef TagIpgLongNormal ipg_type;
  typedef TagNeighborSearchNo neighbor_search_type;
  enum {
    search_type_id = LONG_CUTOFF,
  };
};
struct SEARCH_MODE_GATHER {
  typedef TagForceShort force_type;
  typedef TagSearchShortGather search_type;
  typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
  typedef TagIpgOut ipg_type;
  typedef TagNeighborSearchGather neighbor_search_type;
  enum {
    search_type_id = SHORT_GATHER,
  };
};
struct SEARCH_MODE_SCATTER {
  typedef TagForceShort force_type;
  typedef TagSearchShortScatter search_type;
  typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
  typedef TagIpgIn ipg_type;
  typedef TagNeighborSearchScatter neighbor_search_type;
  enum {
    search_type_id = SHORT_SCATTER,
  };
};
struct SEARCH_MODE_SYMMETRY {
  typedef TagForceShort force_type;
  typedef TagSearchShortSymmetry search_type;
  typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
  typedef TagIpgInAndOut ipg_type;
  typedef TagNeighborSearchSymmetry neighbor_search_type;
  enum {
    search_type_id = SHORT_SYMMETRY,
  };
};

// new TAG for P^3T
struct SEARCH_MODE_LONG_SCATTER {
  typedef TagForceLong force_type;
  typedef TagSearchLongScatter search_type;
  typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
  typedef TagIpgIn ipg_type;
  typedef TagNeighborSearchScatter neighbor_search_type;
  enum {
    search_type_id = LONG_SCATTER,
  };
};
struct SEARCH_MODE_LONG_SYMMETRY {
  typedef TagForceLong force_type;
  typedef TagSearchLongSymmetry search_type;
  typedef TagSearchBoundaryConditionOpenOnly search_boundary_type;
  typedef TagIpgInAndOut ipg_type;
  typedef TagNeighborSearchSymmetry neighbor_search_type;
  enum {
    search_type_id = LONG_SYMMETRY,
  };
};
/*
struct SEARCH_MODE_LONG_CUTOFF_SCATTER{
    typedef TagForceLong force_type;
    typedef TagSearchLongCutoffScatter search_type;
    typedef TagSearchBoundaryConditionOpenPeriodic search_boundary_type;
    typedef TagIpgInAndOut ipg_type;
    enum{
        search_type_id = LONG_CUTOFF_SCATTER,
    };
};
*/

template <class T>
class ValueTypeReduction;
template <>
class ValueTypeReduction<float> {};
template <>
class ValueTypeReduction<double> {};
template <>
class ValueTypeReduction<int> {};
template <>
class ValueTypeReduction<long> {};

enum BOUNDARY_CONDITION {
  BOUNDARY_CONDITION_OPEN,
  BOUNDARY_CONDITION_PERIODIC_X,
  BOUNDARY_CONDITION_PERIODIC_Y,
  BOUNDARY_CONDITION_PERIODIC_Z,
  BOUNDARY_CONDITION_PERIODIC_XY,
  BOUNDARY_CONDITION_PERIODIC_XZ,
  BOUNDARY_CONDITION_PERIODIC_YZ,
  BOUNDARY_CONDITION_PERIODIC_XYZ,
  BOUNDARY_CONDITION_SHEARING_BOX,
  BOUNDARY_CONDITION_USER_DEFINED,
};

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
typedef MPI_Request MpiRequest;
#else
typedef int MpiRequest;
#endif

static const S32 N_SMP_PTCL_TOT_PER_PSYS_DEFAULT = 1000000;

////////
// util
inline F32 CalcSeparationSQPointToBox(const F32vec &point, const F32vec &center,
                                      const F32vec &half_length) {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
  F32 dx = fabs(point.x - center.x);
  F32 dy = fabs(point.y - center.y);
  dx = (half_length.x < dx) ? (dx - half_length.x) : 0.0;
  dy = (half_length.y < dy) ? (dy - half_length.y) : 0.0;
  return dx * dx + dy * dy;
#else
  F32 dx = fabs(point.x - center.x);
  F32 dy = fabs(point.y - center.y);
  F32 dz = fabs(point.z - center.z);
  dx = (half_length.x < dx) ? (dx - half_length.x) : 0.0;
  dy = (half_length.y < dy) ? (dy - half_length.y) : 0.0;
  dz = (half_length.z < dz) ? (dz - half_length.z) : 0.0;
  return dx * dx + dy * dy + dz * dz;
#endif
}

// for check function
inline bool IsInBox(const F32vec &pos, const F32vec &center,
                    const F32 half_length, const F32 tolerance = 1e-6) {
  const F32 tol = -std::abs(tolerance);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
  return ((pos.x - (center.x - half_length)) >= tol) &&
         (((center.x + half_length) - pos.x) >= tol) &&
         ((pos.y - (center.y - half_length)) >= tol) &&
         (((center.y + half_length) - pos.y) >= tol);
#else
  return ((pos.x - (center.x - half_length)) >= tol) &&
         (((center.x + half_length) - pos.x) >= tol) &&
         ((pos.y - (center.y - half_length)) >= tol) &&
         (((center.y + half_length) - pos.y) >= tol) &&
         ((pos.z - (center.z - half_length)) >= tol) &&
         (((center.z + half_length) - pos.z) >= tol);
#endif
}

template <class T>
class Abs {
 public:
  T operator()(const T val) { return std::abs(val); }
};
template <class T>
inline S32 Unique(T val[], const S32 n) {
  S32 ret = 0;
  if (n > 0) {
    ret = 1;
    T ref = val[0];
    for (S32 i = 1; i < n; i++) {
      if (val[i] > ref) {
        ref = val[i];
        val[ret] = ref;
        ret++;
      }
    }
  }
  return ret;
}

template <class T>
inline T GetMSB(const T val);
template <>
inline U64 GetMSB(const U64 val) {
  return (val >> 63) & 0x1;
}
template <>
inline U32 GetMSB(const U32 val) {
  return (val >> 31) & 0x1;
}

template <class T>
inline T ClearMSB(const T val);
template <>
inline U64 ClearMSB(const U64 val) {
  return val & 0x7fffffffffffffff;
}
template <>
inline U32 ClearMSB(const U32 val) {
  return val & 0x7fffffff;
}

template <class T>
inline T SetMSB(const T val);
template <>
inline U64 SetMSB(const U64 val) {
  return val | 0x8000000000000000;
}
template <>
inline U32 SetMSB(const U32 val) {
  return val | 0x80000000;
}

template <class Tp>
inline F64ort GetMinBoxSingleThread(const Tp ptcl[], const S32 n) {
  F64ort pos_box(ptcl[0].getPos(), ptcl[0].getPos());
  for (S32 i = 1; i < n; i++) {
    pos_box.merge(ptcl[i].getPos());
  }
  return pos_box;
}

template <class Tp>
inline F64ort GetMinBox(const Tp ptcl[], const S32 n) {
  F64ort box_loc;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif  // PARTICLE_SIMULATOR_THREAD_PARALLEL
  {
    F64ort box_loc_tmp;
    box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
    for (S32 ip = 0; ip < n; ip++) {
      box_loc_tmp.merge(ptcl[ip].getPos());
    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
    { box_loc.merge(box_loc_tmp); }
  }
  // std::cerr<<"box_loc="<<box_loc<<std::endl;
  const F64vec xlow_loc = box_loc.low_;
  const F64vec xhigh_loc = box_loc.high_;
  const F64vec xlow_glb = Comm::getMinValue(xlow_loc);
  const F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);
  return F64ort(xlow_glb, xhigh_glb);
}

template <class Tp>
inline F64ort GetMinBoxWithMargen(const Tp ptcl[], const S32 n) {
  F64ort box_loc;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
  {
    F64ort box_loc_tmp;
    box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
    for (S32 ip = 0; ip < n; ip++) {
      box_loc_tmp.merge(ptcl[ip].getPos(), ptcl[ip].getRSearch());
    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
    { box_loc.merge(box_loc_tmp); }
  }
  const F64vec xlow_loc = box_loc.low_;
  const F64vec xhigh_loc = box_loc.high_;
  const F64vec xlow_glb = Comm::getMinValue(xlow_loc);
  const F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);
  return F64ort(xlow_glb, xhigh_glb);
}

inline F64 GetWtime() {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
#ifdef PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
  Comm::barrier();
#endif  // PARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
  return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
  // PARTICLE_SIMULATOR_THREAD_PARALLEL
  return omp_get_wtime();
#else
  return (F64)clock() / CLOCKS_PER_SEC;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
}

inline F64 GetWtimeNoBarrier() {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
  return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
  return omp_get_wtime();
#else
  return clock() / CLOCKS_PER_SEC;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
}

struct LessOPForVecX {
  bool operator()(const F64vec &left, const F64vec &right) const {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    return (left.x == right.x)
               ? ((left.y == right.y) ? true : right.y < right.y)
               : left.x < right.x;
#else
    return (left.x == right.x)
               ? ((left.y == right.y)
                      ? ((left.z == right.z) ? true : left.z < right.z)
                      : left.y < right.y)
               : left.x < right.x;
#endif
  }
};

struct LessOPX {
  template <class T>
  bool operator()(const T &left, const T &right) const {
    return left.x < right.x;
  }
};
struct LessOPY {
  template <class T>
  bool operator()(const T &left, const T &right) const {
    return left.y < right.y;
  }
};
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
struct LessOPZ {
  template <class T>
  bool operator()(const T &left, const T &right) const {
    return left.z < right.z;
  }
};
#endif
struct LessOPKEY {
  template <class T>
  bool operator()(const T &left, const T &right) const {
    return left.key_ < right.key_;
  }
};

struct TagRSearch {};
struct TagNoRSearch {};

template <bool T>
struct HasRSearchInner {
  typedef TagRSearch type;
};
template <>
struct HasRSearchInner<false> {
  typedef TagNoRSearch type;
};

template <class T>
class HasRSearch {
 private:
  typedef char One[1];
  typedef char Two[2];

  template <class T2, T2>
  class Check {};

  template <typename T3>
  static One &Func(Check<double (T3::*)(void), &T3::getRSearch> *);
  template <typename T3>
  static One &Func(Check<float (T3::*)(void), &T3::getRSearch> *);
  template <typename T3>
  static One &Func(Check<double (T3::*)(void) const, &T3::getRSearch> *);
  template <typename T3>
  static One &Func(Check<float (T3::*)(void) const, &T3::getRSearch> *);

  template <typename T3>
  static Two &Func(...);

 public:
  typedef typename HasRSearchInner<sizeof(Func<T>(NULL)) == 1>::type type;
};

template <typename Tptcl>
struct HasgetRSearchMethod {
  template <typename U, float (U::*)()>
  struct SFINAE0 {};
  template <typename U, float (U::*)() const>
  struct SFINAE1 {};
  template <typename U, double (U::*)()>
  struct SFINAE2 {};
  template <typename U, double (U::*)() const>
  struct SFINAE3 {};
  template <typename U>
  static char Test(SFINAE0<U, &U::getRSearch> *);
  template <typename U>
  static char Test(SFINAE1<U, &U::getRSearch> *);
  template <typename U>
  static char Test(SFINAE2<U, &U::getRSearch> *);
  template <typename U>
  static char Test(SFINAE3<U, &U::getRSearch> *);
  template <typename U>
  static int Test(...);
  static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
};
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl, std::true_type) {
  return ptcl.getRSearch();
}
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl, std::false_type) {
  return 0.0;
}
template <class Tptcl>
F64 GetMyRSearch(Tptcl ptcl) {
  return GetMyRSearch(
      ptcl, std::integral_constant<bool, HasgetRSearchMethod<Tptcl>::value>());
}

inline F64 GetDistanceMinSq(const F64ort &pos0, const F64ort &pos1,
                            const F64vec &len_peri) {
  const F64vec cen0 = pos0.getCenter();
  const F64vec cen1 = pos1.getCenter();
  const F64vec len = pos0.getHalfLength() + pos1.getHalfLength();
  F64 dis = 0.0;

  F64 dx = fabs(cen0.x - cen1.x);
  dx = (dx < (len_peri.x - dx)) ? dx : (len_peri.x - dx);
  dx = (len.x < dx) ? dx - len.x : 0.0;

  F64 dy = fabs(cen0.y - cen1.y);
  dy = (dy < (len_peri.y - dy)) ? dy : (len_peri.y - dy);
  dy = (len.y < dy) ? dy - len.y : 0.0;
  dis = dx * dx + dy * dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  F64 dz = fabs(cen0.z - cen1.z);
  dz = (dz < (len_peri.z - dz)) ? dz : (len_peri.z - dz);
  dz = (len.z < dz) ? dz - len.z : 0.0;
  dis += dz * dz;
#endif

  return dis;
}

inline F64 GetDistanceMinSq(const F64ort &pos0, const F64vec &pos1,
                            const F64vec &len_peri) {
  const F64vec cen0 = pos0.getCenter();
  const F64vec len = pos0.getHalfLength();
  F64 dis = 0.0;

  F64 dx = fabs(cen0.x - pos1.x);
  dx = (dx < (len_peri.x - dx)) ? dx : (len_peri.x - dx);
  dx = (len.x < dx) ? dx - len.x : 0.0;

  F64 dy = fabs(cen0.y - pos1.y);
  dy = (dy < (len_peri.y - dy)) ? dy : (len_peri.y - dy);
  dy = (len.y < dy) ? dy - len.y : 0.0;

  dis = dx * dx + dy * dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  F64 dz = fabs(cen0.z - pos1.z);
  dz = (dz < (len_peri.z - dz)) ? dz : (len_peri.z - dz);
  dz = (len.z < dz) ? dz - len.z : 0.0;
  dis += dz * dz;
#endif
  return dis;
}

inline F64 GetDistanceMinSq(const F64ort &pos0, const F64ort &pos1) {
  const F64vec cen0 = pos0.getCenter();
  const F64vec cen1 = pos1.getCenter();
  const F64vec len = pos0.getHalfLength() + pos1.getHalfLength();
  F64 dis = 0.0;

  F64 dx = fabs(cen0.x - cen1.x);
  dx = (len.x < dx) ? dx - len.x : 0.0;

  F64 dy = fabs(cen0.y - cen1.y);
  dy = (len.y < dy) ? dy - len.y : 0.0;
  dis = dx * dx + dy * dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  F64 dz = fabs(cen0.z - cen1.z);
  dz = (len.z < dz) ? dz - len.z : 0.0;
  dis += dz * dz;
#endif

  return dis;
}

inline F64 GetDistanceMinSq(const F64ort &pos0, const F64vec &pos1) {
  const F64vec cen0 = pos0.getCenter();
  const F64vec len = pos0.getHalfLength();
  F64 dis = 0.0;

  F64 dx = fabs(cen0.x - pos1.x);
  dx = (len.x < dx) ? dx - len.x : 0.0;

  F64 dy = fabs(cen0.y - pos1.y);
  dy = (len.y < dy) ? dy - len.y : 0.0;

  dis = dx * dx + dy * dy;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  F64 dz = fabs(cen0.z - pos1.z);
  dz = (len.z < dz) ? dz - len.z : 0.0;
  dis += dz * dz;
#endif
  return dis;
}

inline std::string GetBinString(const U32 val) {
  if (!val) return std::string("00000000000000000000000000000000");
  U32 tmp = val;
  std::string str;
  for (S32 i = 0; i < 32; i++) {
    if ((tmp & 1) == 0)
      str.insert(str.begin(), '0');
    else
      str.insert(str.begin(), '1');
    tmp >>= 1;
  }
  return str;
}

inline std::string GetBinString(const U32 *val) {
  if (!*val) return std::string("00000000000000000000000000000000");
  U32 tmp = *val;
  std::string str;
  for (S32 i = 0; i < 32; i++) {
    if ((tmp & 1) == 0)
      str.insert(str.begin(), '0');
    else
      str.insert(str.begin(), '1');
    tmp >>= 1;
  }
  return str;
}

inline std::string GetBinString(const U64 val) {
  if (!val)
    return std::string(
        "0000000000000000000000000000000000000000000000000000000000000000");
  U64 tmp = val;
  std::string str;
  for (S32 i = 0; i < 64; i++) {
    if ((tmp & 1) == 0)
      str.insert(str.begin(), '0');
    else
      str.insert(str.begin(), '1');
    tmp >>= 1;
  }
  return str;
}

inline std::string GetBinString(const U64 *val) {
  if (!*val)
    return std::string(
        "0000000000000000000000000000000000000000000000000000000000000000");
  U64 tmp = *val;
  std::string str;
  for (S32 i = 0; i < 64; i++) {
    if ((tmp & 1) == 0)
      str.insert(str.begin(), '0');
    else
      str.insert(str.begin(), '1');
    tmp >>= 1;
  }
  return str;
}

#ifdef TEST_VARIADIC_TEMPLATE
template <class T>
class IsParticleSystem {
  template <class T2>
  static auto check(T2 x)
      -> decltype(x.ParticleSystemDummyFunc(), std::true_type());
  static auto check(...) -> decltype(std::false_type());

 public:
  typedef decltype(check(std::declval<T>())) value;
};
template <class T>
class IsDomainInfo {
  template <class T2>
  static auto check(T2 x)
      -> decltype(x.DomainInfoDummyFunc(), std::true_type());
  static auto check(...) -> decltype(std::false_type());

 public:
  typedef decltype(check(std::declval<T>())) value;
};
#endif

}  // namespace ParticleSimulator

