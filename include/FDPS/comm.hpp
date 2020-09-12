#ifndef FDPS_COMM_HPP
#define FDPS_COMM_HPP

#include "FDPS/ps_defs.hpp"

namespace ParticleSimulator {

namespace comm_impl {

template <class T>
T allreduceMin(const T &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  T ret;
  T val_tmp = val;
  MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MIN, MPI_COMM_WORLD);
  return ret;
#else
  return val;
#endif
}

template <class T>
T allreduceMax(const T &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  T ret;
  T val_tmp = val;
  MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_MAX, MPI_COMM_WORLD);
  return ret;
#else
  return val;
#endif
}

template <class T>
T allreduceSum(const T &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  T ret;
  T val_tmp = val;
  MPI_Allreduce(&val_tmp, &ret, 1, GetDataType<T>(), MPI_SUM, MPI_COMM_WORLD);
  return ret;
#else
  return val;
#endif
}

template <class T>
void allreduceMin(const T &f_in, const int &i_in, T &f_out,
                         int &i_out) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  struct {
    T x;
    int y;
  } loc, glb;
  loc.x = f_in;
  loc.y = i_in;
  MPI_Allreduce(&loc, &glb, 1, GetDataType<T, int>(), MPI_MINLOC,
                MPI_COMM_WORLD);
  f_out = glb.x;
  i_out = glb.y;
#else
  f_out = f_in;
  i_out = i_in;
#endif
}

template <class T>
void allreduceMax(const T &f_in, const int &i_in, T &f_out,
                         int &i_out) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  struct {
    T x;
    int y;
  } loc, glb;
  loc.x = f_in;
  loc.y = i_in;
  MPI_Allreduce(&loc, &glb, 1, GetDataType<T, int>(), MPI_MAXLOC,
                MPI_COMM_WORLD);
  f_out = glb.x;
  i_out = glb.y;
#else
  f_out = f_in;
  i_out = i_in;
#endif
}

template <typename T>
struct minValue {
  static T get(const T &val) { return allreduceMin(val); }
};

template <>
struct minValue<Vector2<float>> {
  static Vector2<float> get(const Vector2<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<float> ret;
    MPI_Allreduce((float *)&val.x, (float *)&ret.x, 2, MPI_FLOAT, MPI_MIN,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct minValue<Vector2<double>> {
  static Vector2<double> get(const Vector2<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<double> ret;
    MPI_Allreduce((double *)&val.x, (double *)&ret.x, 2, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct minValue<Vector3<float>> {
  static Vector3<float> get(const Vector3<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<float> ret;
    MPI_Allreduce((float *)&val.x, (float *)&ret.x, 3, MPI_FLOAT, MPI_MIN,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct minValue<Vector3<double>> {
  static Vector3<double> get(const Vector3<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<double> ret;
    MPI_Allreduce((double *)&val.x, (double *)&ret.x, 3, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <typename T>
struct maxValue {
  static T get(const T &val) { return allreduceMax(val); }
};

template <>
struct maxValue<Vector2<float>> {
  static Vector2<float> get(const Vector2<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<float> ret;
    MPI_Allreduce((float *)&val.x, (float *)&ret.x, 2, MPI_FLOAT, MPI_MAX,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct maxValue<Vector2<double>> {
  static Vector2<double> get(const Vector2<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<double> ret;
    MPI_Allreduce((double *)&val.x, (double *)&ret.x, 2, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct maxValue<Vector3<float>> {
  static Vector3<float> get(const Vector3<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<float> ret;
    MPI_Allreduce((float *)&val.x, (float *)&ret.x, 3, MPI_FLOAT, MPI_MAX,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct maxValue<Vector3<double>> {
  static Vector3<double> get(const Vector3<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<double> ret;
    MPI_Allreduce((double *)&val.x, (double *)&ret.x, 3, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <typename T>
struct sum {
  static T get(const T &val) { return allreduceSum(val); }
};

template <>
struct sum<Vector3<float>> {
  static Vector3<float> get(const Vector3<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<float> ret;
    MPI_Allreduce((float *)&val, (float *)&ret, 3, MPI_FLOAT, MPI_SUM,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct sum<Vector3<double>> {
  static Vector3<double> get(const Vector3<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector3<double> ret;
    MPI_Allreduce((double *)&val, (double *)&ret, 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct sum<Vector2<float>> {
  static Vector2<float> get(const Vector2<float> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<float> ret;
    MPI_Allreduce((float *)&val, (float *)&ret, 2, MPI_FLOAT, MPI_SUM,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

template <>
struct sum<Vector2<double>> {
  static Vector2<double> get(const Vector2<double> &val) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    Vector2<double> ret;
    MPI_Allreduce((double *)&val, (double *)&ret, 2, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return ret;
#else
    return val;
#endif
  }
};

}  // namespace comm_impl

class Comm {
 private:
  S32 rank_;
  S32 n_proc_;
  S32 n_thread_;
  S32 rank_multi_dim_[DIMENSION];
  S32 n_proc_multi_dim_[DIMENSION];

  Comm();

  static Comm &getInstance();

 public:
  static void setNumberOfProcMultiDim(const S32 id, const S32 n) {
    getInstance().n_proc_multi_dim_[id] = n;
  }

  static void setRankMultiDim(const S32 id, const S32 r) {
    getInstance().rank_multi_dim_[id] = r;
  }

  static S32 getRank() { return Comm::getInstance().rank_; }

  static S32 getNumberOfProc() { return getInstance().n_proc_; }

  static S32 getRankMultiDim(const S32 id) {
    return getInstance().rank_multi_dim_[id];
  }

  static S32 getNumberOfProcMultiDim(const S32 id) {
    return getInstance().n_proc_multi_dim_[id];
  }

  static S32 getNumberOfThread() { return getInstance().n_thread_; }

  static S32 getNumThreads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
    return omp_get_num_threads();
#else
    return 1;
#endif
  }

  static S32 getThreadNum() {
    //#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
    return omp_get_thread_num();
#else
    return 0;
#endif
  }

  static void barrier() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  static bool synchronizeConditionalBranchAND(const bool &local) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    bool global;
    bool local_tmp = local;
    MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    return global;
#else
    return local;
#endif
  }

  static bool synchronizeConditionalBranchOR(const bool &local) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    bool global;
    bool local_tmp = local;
    // std::cerr<<"rank="<<Comm::getRank()<<" local="<<local<<std::endl;
    MPI_Allreduce(&local_tmp, &global, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    // std::cerr<<"rank="<<Comm::getRank()<<" global="<<global<<std::endl;
    return global;
#else
    return local;
#endif
  }

  ///////////////////////////
  // MPI ALLREDUCE WRAPPER //
  template <class T>
  static T getMinValue(const T &val) {
    return comm_impl::minValue<T>::get(val);
  }

  template <class T>
  static T getMaxValue(const T &val) {
    return comm_impl::maxValue<T>::get(val);
  }

  template <class Tfloat>
  static void getMinValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out,
                          int &i_out) {
    comm_impl::allreduceMin(f_in, i_in, f_out, i_out);
  }

  template <class Tfloat>
  static void getMaxValue(const Tfloat &f_in, const int &i_in, Tfloat &f_out,
                          int &i_out) {
    comm_impl::allreduceMax(f_in, i_in, f_out, i_out);
  }

  template <class T>
  static T getSum(const T &val) {
    return comm_impl::sum<T>::get(val);
  }

  ///////////////////////
  // MPI BCAST WRAPPER //
  // new functions 10 Feb 2015
  template <class T>
  static inline void broadcast(T *val, const int n, const int src = 0) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Bcast(val, n, GetDataType<T>(), src, MPI_COMM_WORLD);
#else
    // NOP
#endif
  }

  ///////////////////////////
  // MPI GATHER WRAPPER //
  template <class T>
  static inline void gather(T *val_send, int n, T *val_recv) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Gather(val_send, n, GetDataType<T>(), val_recv, n, GetDataType<T>(), 0,
               MPI_COMM_WORLD);
#else
    for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
  }

  template <class T>
  static inline void gatherV(T *val_send, int n_send, T *val_recv) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    const int n_proc = Comm::getNumberOfProc();
    int *n_recv = new int[n_proc];
    int *n_disp_recv = new int[n_proc + 1];
    gather(&n_send, 1, n_recv);
    n_disp_recv[0] = 0;
    for (S32 i = 0; i < n_proc; i++) {
      n_disp_recv[i + 1] = n_disp_recv[i] + n_recv[i];
    }
    MPI_Gatherv(val_send, n_send, GetDataType<T>(), val_recv, n_recv,
                n_disp_recv, GetDataType<T>(), 0, MPI_COMM_WORLD);
    delete[] n_recv;
    delete[] n_disp_recv;
#else
    int n = n_send;
    for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
  }

  template <class T>
  static inline void gatherV(T *val_send,      // in
                             int n_send,       // in
                             T *val_recv,      // out
                             int *n_recv,      // in
                             int *n_recv_disp  // in
  ) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Gatherv(val_send, n_send, GetDataType<T>(), val_recv, n_recv,
                n_recv_disp, GetDataType<T>(), 0, MPI_COMM_WORLD);
#else
    for (int i = 0; i < n_send; i++) val_recv[i] = val_send[i];
      // n_recv[0] = n_recv_disp[1] = n_send; //not needed ?
      // n_recv_disp[0] = 0; //not needed ?
#endif
  }

  ///////////////////////////
  // MPI SCATTER WRAPPER //
  template <class T>
  static inline void scatterV(T *val_send, int *n_send, int *n_send_disp,
                              T *val_recv,  // output
                              int n_recv, int root = 0) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Scatterv(val_send, n_send, n_send_disp, GetDataType<T>(), val_recv,
                 n_recv, GetDataType<T>(), root, MPI_COMM_WORLD);
#else
    S32 n_proc = Comm::getNumberOfProc();
    for (int i = 0; i < n_send_disp[n_proc]; i++) val_recv[i] = val_send[i];
#endif
  }

  ///////////////////////////
  // MPI ALLGATHER WRAPPER //
  template <class T>
  static inline void allGather(const T *val_send,  // in
                               const int n,        // in
                               T *val_recv         // out
  ) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    T *val_send_tmp = const_cast<T *>(val_send);
    MPI_Allgather(val_send_tmp, n, GetDataType<T>(), val_recv, n,
                  GetDataType<T>(), MPI_COMM_WORLD);
#else
    for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
  }

  template <class T>
  static inline void allGatherV(const T *val_send,  // in
                                const int n_send,   // in
                                T *val_recv,        // out
                                int *n_recv,        // in
                                int *n_recv_disp    // in
  ) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    T *val_send_tmp = const_cast<T *>(val_send);
    int *n_recv_tmp = const_cast<int *>(n_recv);
    int *n_recv_disp_tmp = const_cast<int *>(n_recv_disp);
    MPI_Allgatherv(val_send_tmp, n_send, GetDataType<T>(), val_recv, n_recv_tmp,
                   n_recv_disp_tmp, GetDataType<T>(), MPI_COMM_WORLD);
#else
    for (int i = 0; i < n_send; i++) val_recv[i] = val_send[i];
    n_recv[0] = n_recv_disp[1] = n_send;
    n_recv_disp[0] = 0;
#endif
  }

  ///////////////////////////
  // MPI ALLTOALL WRAPPER //
  template <class T>
  static inline void allToAll(const T *val_send, const int n, T *val_recv) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    T *val_send_tmp = const_cast<T *>(val_send);
    MPI_Alltoall(val_send_tmp, n, GetDataType<T>(), val_recv, n,
                 GetDataType<T>(), MPI_COMM_WORLD);
#else
    for (int i = 0; i < n; i++) val_recv[i] = val_send[i];
#endif
  }

  template <class T>
  static inline void allToAllV(const T *val_send, const int *n_send,
                               const int *n_send_disp, T *val_recv, int *n_recv,
                               int *n_recv_disp) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    T *val_send_tmp = const_cast<T *>(val_send);
    int *n_send_tmp = const_cast<int *>(n_send);
    int *n_send_disp_tmp = const_cast<int *>(n_send_disp);
    int *n_recv_tmp = const_cast<int *>(n_recv);
    int *n_recv_disp_tmp = const_cast<int *>(n_recv_disp);
    MPI_Alltoallv(val_send_tmp, n_send_tmp, n_send_disp_tmp, GetDataType<T>(),
                  val_recv, n_recv_tmp, n_recv_disp_tmp, GetDataType<T>(),
                  MPI_COMM_WORLD);
#else
    for (int i = 0; i < n_send[0]; i++) val_recv[i] = val_send[i];
    n_recv[0] = n_send[0];
    n_recv_disp[0] = n_send_disp[0];
    n_recv_disp[1] = n_send_disp[1];
#endif
  }

};  // END OF Comm

}  // namespace ParticleSimulator

#endif  // FDPS_COMM_HPP