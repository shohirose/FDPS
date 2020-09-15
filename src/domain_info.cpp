#include "FDPS/domain_info.hpp"

namespace ParticleSimulator {

void DomainInfo::calculateBoundaryOfDomainX(const S32 np,
                                            const F64vec pos_sample[],
                                            const S32 istart, const S32 iend,
                                            F64 &xlow, F64 &xhigh) {
  if (istart == 0)
    xlow = pos_root_domain_.low_.x;
  else
    xlow = 0.5 * (pos_sample[istart - 1].x + pos_sample[istart].x);
  if (iend == np - 1)
    xhigh = pos_root_domain_.high_.x;
  else
    xhigh = 0.5 * (pos_sample[iend].x + pos_sample[iend + 1].x);
}

void DomainInfo::calculateBoundaryOfDomainY(const S32 np,
                                            const F64vec pos_sample[],
                                            const S32 istart, const S32 iend,
                                            F64 &xlow, F64 &xhigh) {
  if (istart == 0)
    xlow = pos_root_domain_.low_.y;
  else
    xlow = 0.5 * (pos_sample[istart - 1].y + pos_sample[istart].y);
  if (iend == np - 1)
    xhigh = pos_root_domain_.high_.y;
  else
    xhigh = 0.5 * (pos_sample[iend].y + pos_sample[iend + 1].y);
}

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
void DomainInfo::calculateBoundaryOfDomainZ(const S32 np,
                                            const F64vec pos_sample[],
                                            const S32 istart, const S32 iend,
                                            F64 &xlow, F64 &xhigh) {
  if (istart == 0)
    xlow = pos_root_domain_.low_.z;
  else
    xlow = 0.5 * (pos_sample[istart - 1].z + pos_sample[istart].z);
  if (iend == np - 1)
    xhigh = pos_root_domain_.high_.z;
  else
    xhigh = 0.5 * (pos_sample[iend].z + pos_sample[iend + 1].z);
}
#endif

DomainInfo::DomainInfo() {
  first_call_by_initialize = true;
  first_call_by_decomposeDomain = true;
  periodic_axis_[0] = periodic_axis_[1] = false;
  pos_root_domain_.low_.x = -LARGE_FLOAT;
  pos_root_domain_.high_.x = LARGE_FLOAT;
  pos_root_domain_.low_.y = -LARGE_FLOAT;
  pos_root_domain_.high_.y = LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  periodic_axis_[2] = false;
  pos_root_domain_.low_.z = -LARGE_FLOAT;
  pos_root_domain_.high_.z = LARGE_FLOAT;
#endif
  boundary_condition_ = BOUNDARY_CONDITION_OPEN;
}

DomainInfo::~DomainInfo() {
  delete[] n_smp_array_;
  delete[] n_smp_disp_array_;
  delete[] pos_domain_;
  delete[] pos_domain_temp_;
  delete[] pos_sample_tot_;
  delete[] pos_sample_loc_;
}

void DomainInfo::initialize(const F32 coef_ema) {
  if (coef_ema < 0.0 || coef_ema > 1.0) {
    PARTICLE_SIMULATOR_PRINT_ERROR(
        "The smoothing factor of an exponential moving average is must "
        "between 0 and 1.");
    std::cerr << "The smoothing factor of an exponential moving average is "
                 "must between 0 and 1."
              << std::endl;
    Abort(-1);
  }
  assert(first_call_by_initialize);

  first_call_by_initialize = false;
  pos_sample_tot_ = NULL;
  pos_sample_loc_ = NULL;

  n_smp_array_ = new S32[Comm::getNumberOfProc()];

  n_smp_disp_array_ = new S32[Comm::getNumberOfProc() + 1];

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  pos_domain_ = new F64ort[Comm::getNumberOfProc()];
  pos_domain_temp_ = new F64ort[Comm::getNumberOfProc()];
#else
  pos_domain_ = new F64ort[1];
  pos_domain_temp_ = new F64ort[1];
#endif

  coef_ema_ = coef_ema;
  target_number_of_sample_particle_ = 0;
  number_of_sample_particle_tot_ = 0;
  number_of_sample_particle_loc_ = 0;

  // S32 rank_tmp[DIMENSION];
  S32 rank_tmp[DIMENSION_LIMIT];
  SetNumberOfDomainMultiDimension<DIMENSION>(n_domain_, rank_tmp);

  // std::cerr<<"check 2"<<std::endl;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  // NEW
  int rank_glb = Comm::getRank();

  for (S32 d = DIMENSION - 1; d >= 0; d--) {
    rank_1d_[d] = rank_glb % n_domain_[d];
    rank_glb /= n_domain_[d];
    MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_ + d);
    MPI_Comm_rank(comm_sub_[d], rank_sub_ + d);
    MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_ + d);
    MPI_Comm_size(comm_sub_[d], n_proc_sub_ + d);
  }

  for (S32 d = DIMENSION - 1; d >= 0; d--) {
    Comm::setRankMultiDim(d, rank_tmp[d]);
    Comm::setNumberOfProcMultiDim(d, n_domain_[d]);
  }

#endif
}

void DomainInfo::setNumberOfDomainMultiDimension(const S32 nx, const S32 ny,
                                                 const S32 nz) {
  S32 n_proc = Comm::getNumberOfProc();
  if (n_proc != nx * ny * nz) {
    PARTICLE_SIMULATOR_PRINT_ERROR(
        "devided number of domains is not consistent with total processe "
        "number");
    Abort(-1);
  }
  n_domain_[0] = nx;
  n_domain_[1] = ny;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
  n_domain_[2] = 1;
#else
  n_domain_[2] = nz;
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  int rank_glb = Comm::getRank();
  // make comunicator
  for (S32 d = DIMENSION - 1; d >= 0; d--) {
    rank_1d_[d] = rank_glb % n_domain_[d];
    rank_glb /= n_domain_[d];
    MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_ + d);
    MPI_Comm_rank(comm_sub_[d], rank_sub_ + d);
    MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_ + d);
    MPI_Comm_size(comm_sub_[d], n_proc_sub_ + d);
  }
  int rank_tmp = Comm::getRank();
  for (S32 d = DIMENSION - 1; d >= 0; d--) {
    Comm::setRankMultiDim(d, rank_tmp % n_domain_[d]);
    rank_tmp /= n_domain_[d];
  }
  Comm::setNumberOfProcMultiDim(0, nx);
  Comm::setNumberOfProcMultiDim(1, ny);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  Comm::setNumberOfProcMultiDim(2, ny);
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
}

// new version multi-dimensional gathering
void DomainInfo::decomposeDomainMultiStep() {
  F64 time_offset = GetWtime();
  // assert(!first_call_by_decomposeDomain);
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL  // ifNdef
  pos_domain_[0] = pos_root_domain_;
#else
  static bool first = true;
  static S32 *n_send;
  static S32 *n_recv;
  static S32 *n_send_disp;
  static S32 *n_recv_disp;
  static S32 *i_head;
  static S32 *i_tail;
  static F64vec *pos_sample_buf;
  static F64 *coord_buf;
  static F64 *coord_tot;
  static F64 *x_coord;
  static F64 *y_coord;
  static F64ort *pos_domain_temp_buf;
  S32 n_proc_glb = Comm::getNumberOfProc();
  // S32 rank_glb = Comm::getRank();

  if (first) {
    n_send = new S32[n_proc_glb];
    n_recv = new S32[n_proc_glb];
    n_send_disp = new S32[n_proc_glb + 1];
    n_recv_disp = new S32[n_proc_glb + 1];
    i_head = new S32[n_proc_glb];
    i_tail = new S32[n_proc_glb];
    pos_sample_buf = new F64vec[target_number_of_sample_particle_];
    coord_buf = new F64[n_proc_glb * 2];
    coord_tot = new F64[n_proc_glb * 2];
    x_coord = new F64[n_proc_glb + 1];
    y_coord = new F64[n_proc_glb + 1];
    pos_domain_temp_buf = new F64ort[n_proc_glb];
    first = false;
  }

  ///////////// sort particles along x direction
  std::sort(pos_sample_loc_, pos_sample_loc_ + number_of_sample_particle_loc_,
            LessOPX());

  ///////////// migrate particles along x direction
  for (S32 i = 0; i < n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
  S32 id_domain_3d = 0;
  S32 id_domain_x = 0;
  for (S32 i = 0; i < number_of_sample_particle_loc_; i++) {
    while (pos_domain_[id_domain_3d].high_.x <= pos_sample_loc_[i].x) {
      id_domain_3d += n_proc_sub_[0];
      id_domain_x++;
    }
    n_send[id_domain_x]++;
  }
  MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(),
               comm_1d_[0]);
  n_send_disp[0] = n_recv_disp[0] = 0;
  for (S32 i = 0; i < n_domain_[0]; i++) {
    n_send_disp[i + 1] = n_send_disp[i] + n_send[i];
    n_recv_disp[i + 1] = n_recv_disp[i] + n_recv[i];
  }
  MPI_Alltoallv(pos_sample_loc_, n_send, n_send_disp, GetDataType<F64vec>(),
                pos_sample_buf, n_recv, n_recv_disp, GetDataType<F64vec>(),
                comm_1d_[0]);

  ///////////// allgather particles in Y-Z plane
  S32 n_send_tmp = n_recv_disp[n_domain_[0]];  // # of particle in own cell.
  MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1,
                GetDataType<S32>(), comm_sub_[0]);
  n_recv_disp[0] = 0;
  for (S32 i = 0; i < n_proc_sub_[0]; i++) {
    n_recv_disp[i + 1] = n_recv_disp[i] + n_recv[i];
  }
  S32 n_par_slab = n_recv_disp[n_proc_sub_[0]];

  MPI_Allgatherv(pos_sample_buf, n_send_tmp, GetDataType<F64vec>(),
                 pos_sample_tot_, n_recv, n_recv_disp, GetDataType<F64vec>(),
                 comm_sub_[0]);

  ///////////// sort particles along x direction again
  std::sort(pos_sample_tot_, pos_sample_tot_ + n_par_slab, LessOPX());
  ///////////// determine X coord
  MPI_Allgather(&n_par_slab, 1, GetDataType<S32>(), n_recv, 1,
                GetDataType<S32>(), comm_1d_[0]);
  n_recv_disp[0] = 0;
  for (S32 i = 0; i < n_domain_[0]; i++) {
    n_recv_disp[i + 1] = n_recv_disp[i] + n_recv[i];
  }
  number_of_sample_particle_tot_ = n_recv_disp[n_domain_[0]];

  // get index of
  S32 n_ave = number_of_sample_particle_tot_ / n_domain_[0];
  for (S32 i = 0; i < n_domain_[0]; i++) {
    i_head[i] = n_ave * i;
    if (i < number_of_sample_particle_tot_ % n_domain_[0]) {
      i_head[i] += i;
    } else {
      i_head[i] += number_of_sample_particle_tot_ % n_domain_[0];
    }
    if (i > 0) i_tail[i - 1] = i_head[i] - 1;
  }
  i_tail[n_domain_[0] - 1] = number_of_sample_particle_tot_ - 1;

  n_send_tmp = 0;  // temporally used
  for (S32 i = 0; i < n_domain_[0]; i++) {
    if (n_recv_disp[rank_1d_[0]] <= i_head[i] &&
        i_head[i] < n_recv_disp[rank_1d_[0]] + n_par_slab) {
      S32 i_tmp = i_head[i] - n_recv_disp[rank_1d_[0]];
      coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
    }
    if (n_recv_disp[rank_1d_[0]] <= i_tail[i] &&
        i_tail[i] < n_recv_disp[rank_1d_[0]] + n_par_slab) {
      S32 i_tmp = i_tail[i] - n_recv_disp[rank_1d_[0]];
      coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
    }
  }

  MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1,
                GetDataType<S32>(), comm_1d_[0]);
  n_recv_disp[0] = 0;
  for (S32 i = 0; i < n_domain_[0]; i++) {
    n_recv_disp[i + 1] = n_recv_disp[i] + n_recv[i];
  }

  MPI_Allgatherv(coord_buf, n_send_tmp, GetDataType<>(coord_buf[0]), coord_tot,
                 n_recv, n_recv_disp, GetDataType<>(coord_buf[0]), comm_1d_[0]);

  assert(n_recv_disp[n_domain_[0]] == n_domain_[0] * 2);

  // size of x_coord_buf is n_domain_[0]+1
  x_coord[0] = pos_root_domain_.low_.x;
  x_coord[n_domain_[0]] = pos_root_domain_.high_.x;

  for (S32 i = 1; i < n_domain_[0]; i++) {
    x_coord[i] = (coord_tot[i * 2] + coord_tot[i * 2 - 1]) * 0.5;
  }

  ///////////// migrate particles along x direction
  for (S32 i = 0; i < n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
  id_domain_x = 0;
  for (S32 i = 0; i < n_par_slab; i++) {
    while (x_coord[id_domain_x + 1] <= pos_sample_tot_[i].x) id_domain_x++;
    n_send[id_domain_x]++;
  }

  MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(),
               comm_1d_[0]);
  n_send_disp[0] = n_recv_disp[0] = 0;
  for (S32 i = 0; i < n_domain_[0]; i++) {
    n_send_disp[i + 1] = n_send_disp[i] + n_send[i];
    n_recv_disp[i + 1] = n_recv_disp[i] + n_recv[i];
  }
  MPI_Alltoallv(pos_sample_tot_, n_send, n_send_disp, GetDataType<F64vec>(),
                pos_sample_buf, n_recv, n_recv_disp, GetDataType<F64vec>(),
                comm_1d_[0]);
  n_par_slab = n_recv_disp[n_domain_[0]];

  ////////////////////////////////////
  ///////////// determine y corrdinate
  std::sort(pos_sample_buf, pos_sample_buf + n_par_slab, LessOPY());

  // get index of
  n_ave = n_par_slab / n_domain_[1];
  for (S32 i = 0; i < n_domain_[1]; i++) {
    i_head[i] = n_ave * i;
    if (i < n_par_slab % n_domain_[1])
      i_head[i] += i;
    else
      i_head[i] += n_par_slab % n_domain_[1];
    if (i > 0) i_tail[i - 1] = i_head[i] - 1;
  }
  i_tail[n_domain_[1] - 1] = n_par_slab - 1;

  // size of y_coord is n_domain_[1]+1
  y_coord[0] = pos_root_domain_.low_.y;
  y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
  for (S32 i = 1; i < n_domain_[1]; i++)
    y_coord[i] =
        (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i - 1]].y) * 0.5;

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  ////////////////////////////////////
  ///////////// determine z corrdinate
  //#pragma omp parallel for
  for (S32 iy = 0; iy < n_domain_[1]; iy++) {
    const S32 iy_ptcl_head = i_head[iy];
    const S32 iy_ptcl_tail = i_tail[iy];
    const S32 nz_tot = iy_ptcl_tail - iy_ptcl_head + 1;
    std::sort(pos_sample_buf + iy_ptcl_head,
              pos_sample_buf + iy_ptcl_head + nz_tot, LessOPZ());

    S32 nz_ave_tmp = nz_tot / n_domain_[2];
    pos_domain_temp_buf[iy * n_domain_[2]].low_.z = pos_root_domain_.low_.z;
    pos_domain_temp_buf[(iy + 1) * n_domain_[2] - 1].high_.z =
        pos_root_domain_.high_.z;
    pos_domain_temp_buf[iy * n_domain_[2]].low_.y = y_coord[iy];
    pos_domain_temp_buf[iy * n_domain_[2]].high_.y = y_coord[iy + 1];
    for (S32 iz = 1; iz < n_domain_[2]; iz++) {
      pos_domain_temp_buf[iy * n_domain_[2] + iz].low_.y = y_coord[iy];
      pos_domain_temp_buf[iy * n_domain_[2] + iz].high_.y = y_coord[iy + 1];
      S32 iz_tmp = nz_ave_tmp * iz;
      if (iz < nz_tot % n_domain_[2])
        iz_tmp += iz;
      else
        iz_tmp += nz_tot % n_domain_[2];
      F64 z_coord_tmp = (pos_sample_buf[iy_ptcl_head + iz_tmp].z +
                         pos_sample_buf[iy_ptcl_head + iz_tmp - 1].z) *
                        0.5;

      pos_domain_temp_buf[iy * n_domain_[2] + iz].low_.z = z_coord_tmp;
      pos_domain_temp_buf[iy * n_domain_[2] + iz - 1].high_.z = z_coord_tmp;
    }
  }
#endif
  for (S32 i = 0; i < n_proc_sub_[0]; i++) {
    pos_domain_temp_buf[i].low_.x = x_coord[rank_1d_[0]];
    pos_domain_temp_buf[i].high_.x = x_coord[rank_1d_[0] + 1];
  }

  //////////////////////////////////////////////
  ///////////// exchange pos_domain_tmp
  MPI_Allgather(pos_domain_temp_buf, n_proc_sub_[0], GetDataType<F64ort>(),
                pos_domain_temp_, n_proc_sub_[0], GetDataType<F64ort>(),
                comm_1d_[0]);

  if (first_call_by_decomposeDomain) {
    first_call_by_decomposeDomain = false;
    for (S32 i = 0; i < n_proc_glb; i++) {
      pos_domain_[i].low_ = pos_domain_temp_[i].low_;
      pos_domain_[i].high_ = pos_domain_temp_[i].high_;
    }
  } else {
    for (S32 i = 0; i < n_proc_glb; i++) {
      pos_domain_[i].low_ = (F64)coef_ema_ * pos_domain_temp_[i].low_ +
                            (F64)(1. - coef_ema_) * pos_domain_[i].low_;
      pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ +
                             (F64)(1. - coef_ema_) * pos_domain_[i].high_;
    }
  }
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
  time_profile_.decompose_domain = GetWtime() - time_offset;
}

void DomainInfo::decomposeDomain() {
  F64 time_offset = GetWtime();
  // ****** collect sample particles to process 0. ******
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  S32 nproc = Comm::getNumberOfProc();
  S32 myrank = Comm::getRank();
#ifdef __HPC_ACE__
  Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
  n_smp_disp_array_[0] = 0;
  for (S32 i = 0; i < nproc; i++) {
    n_smp_disp_array_[i + 1] = n_smp_disp_array_[i] + n_smp_array_[i];
  }
  Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_,
                   pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
  number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#else   //__HPC_ACE__
  Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);

  n_smp_disp_array_[0] = 0;
  for (S32 i = 0; i < nproc; i++) {
    n_smp_disp_array_[i + 1] = n_smp_disp_array_[i] + n_smp_array_[i];
  }

  Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_,
                   pos_sample_tot_, n_smp_array_, n_smp_disp_array_);

  number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#endif  //__HPC_ACE__

  // ****************************************************
  // *** decompose domain *******************************
  if (myrank == 0) {
    S32 *istart = new S32[nproc];
    S32 *iend = new S32[nproc];
    // --- x direction --------------------------
    std::sort(pos_sample_tot_, pos_sample_tot_ + number_of_sample_particle_tot_,
              Cmpvec(&F64vec::x));
    for (S32 i = 0; i < nproc; i++) {
      istart[i] =
          ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
      if (i > 0) iend[i - 1] = istart[i] - 1;
    }
    iend[nproc - 1] = number_of_sample_particle_tot_ - 1;
    for (S32 ix = 0; ix < n_domain_[0]; ix++) {
      S32 ix0 = ix * n_domain_[1] * n_domain_[2];
      S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
      F64 x0 = 0.0;
      F64 x1 = 0.0;
      calculateBoundaryOfDomainX(number_of_sample_particle_tot_,
                                 pos_sample_tot_, istart[ix0], iend[ix1 - 1],
                                 x0, x1);
      for (S32 i = ix0; i < ix1; i++) {
        pos_domain_temp_[i].low_.x = x0;
        pos_domain_temp_[i].high_.x = x1;
      }
    }
    // ------------------------------------------
    // --- y direction --------------------------
    for (S32 ix = 0; ix < n_domain_[0]; ix++) {
      S32 ix0 = ix * n_domain_[1] * n_domain_[2];
      S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
      std::sort(pos_sample_tot_ + istart[ix0],
                pos_sample_tot_ + (iend[ix1 - 1] + 1), Cmpvec(&F64vec::y));
      S32 number_of_sample_particle_tot_y = iend[ix1 - 1] - istart[ix0] + 1;
      for (S32 iy = 0; iy < n_domain_[1]; iy++) {
        S32 iy0 = ix0 + iy * n_domain_[2];
        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
        // F64 y0, y1;
        F64 y0 = 0.0;
        F64 y1 = 0.0;

        calculateBoundaryOfDomainY(
            number_of_sample_particle_tot_y, pos_sample_tot_ + istart[ix0],
            istart[iy0] - istart[ix0], iend[iy1 - 1] - istart[ix0], y0, y1);
        for (S32 i = iy0; i < iy1; i++) {
          pos_domain_temp_[i].low_.y = y0;
          pos_domain_temp_[i].high_.y = y1;
        }
      }
    }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    // ------------------------------------------
    // --- z direction --------------------------
    for (S32 ix = 0; ix < n_domain_[0]; ix++) {
      S32 ix0 = ix * n_domain_[1] * n_domain_[2];
      for (S32 iy = 0; iy < n_domain_[1]; iy++) {
        S32 iy0 = ix0 + iy * n_domain_[2];
        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];

        std::sort(pos_sample_tot_ + istart[iy0],
                  pos_sample_tot_ + (iend[iy1 - 1] + 1), Cmpvec(&F64vec::z));
        S32 number_of_sample_particle_tot_z = iend[iy1 - 1] - istart[iy0] + 1;
        for (S32 iz = 0; iz < n_domain_[2]; iz++) {
          S32 iz0 = iy0 + iz;
          // F64 z0, z1;
          F64 z0 = 0.0;
          F64 z1 = 0.0;

          calculateBoundaryOfDomainZ(
              number_of_sample_particle_tot_z, pos_sample_tot_ + istart[iy0],
              istart[iz0] - istart[iy0], iend[iz0] - istart[iy0], z0, z1);
          pos_domain_temp_[iz0].low_.z = z0;
          pos_domain_temp_[iz0].high_.z = z1;
        }
      }
    }
#endif  // PARTICLE_SIMULATOR_TWO_DIMENSION
        // ------------------------------------------
        // --- process first ------------------------
    if (first_call_by_decomposeDomain) {
      for (S32 i = 0; i < nproc; i++) {
        pos_domain_[i].low_ = pos_domain_temp_[i].low_;
        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
      }
    } else {
      for (S32 i = 0; i < nproc; i++) {
        pos_domain_[i].low_ = (F64)coef_ema_ * pos_domain_temp_[i].low_ +
                              (F64)(1. - coef_ema_) * pos_domain_[i].low_;
        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ +
                               (F64)(1. - coef_ema_) * pos_domain_[i].high_;
      }
    }
    // ------------------------------------------
    delete[] istart;
    delete[] iend;
  }
  // ****************************************************
  // *** broad cast pos_domain_ *************************
  MPI_Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0, MPI_COMM_WORLD);
  if (first_call_by_decomposeDomain) {
    first_call_by_decomposeDomain = false;
    MPI_Bcast(&first_call_by_decomposeDomain, 1, GetDataType<bool>(), 0,
              MPI_COMM_WORLD);
  }
  // ****************************************************
#else   // PARTICLE_SIMULATOR_MPI_PARALLEL
  pos_domain_[0] = pos_root_domain_;
#endif  // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
  PARTICLE_SIMULATOR_PRINT_LINE_INFO();
  std::cout << "pos_root_domain_=" << pos_root_domain_ << std::endl;
  std::cout << "pos_domain_[Comm::getRank()]=" << pos_domain_[Comm::getRank()]
            << std::endl;
#endif
  time_profile_.decompose_domain += GetWtime() - time_offset;
}

void DomainInfo::setBoundaryCondition(enum BOUNDARY_CONDITION bc) {
  boundary_condition_ = bc;
  if (DIMENSION == 2 && (bc == BOUNDARY_CONDITION_PERIODIC_XYZ ||
                         bc == BOUNDARY_CONDITION_PERIODIC_XZ ||
                         bc == BOUNDARY_CONDITION_PERIODIC_YZ ||
                         bc == BOUNDARY_CONDITION_PERIODIC_Z)) {
    throw "PS_ERROR: in setBoundaryCondition(enum BOUNDARY_CONDITION) \n boundary condition is incompatible with DIMENSION";
  }
  if (bc == BOUNDARY_CONDITION_PERIODIC_X)
    periodic_axis_[0] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_Y)
    periodic_axis_[1] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_Z)
    periodic_axis_[2] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_XY)
    periodic_axis_[0] = periodic_axis_[1] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_XZ)
    periodic_axis_[0] = periodic_axis_[2] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_YZ)
    periodic_axis_[1] = periodic_axis_[2] = true;
  else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ)
    periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = true;
}

void DomainInfo::setPosRootDomain(const F64vec &low, const F64vec &high) {
  for (S32 k = 0; k < DIMENSION; k++) {
    if (low[k] > high[k]) {
      PARTICLE_SIMULATOR_PRINT_ERROR(
          "The coodinate of the root domain is inconsistent.");
      std::cerr << "The coordinate of the low vertex of the rood domain=" << low
                << std::endl;
      std::cerr << "The coordinate of the high vertex of the rood domain="
                << high << std::endl;
      Abort(-1);
    }
  }

  for (S32 i = 0; i < DIMENSION; i++) {
    if (periodic_axis_[i] == false) continue;
    if (low[i] < high[i]) {
      pos_root_domain_.low_[i] = low[i];
      pos_root_domain_.high_[i] = high[i];
    }
  }
}

}  // namespace ParticleSimulator