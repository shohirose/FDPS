#pragma once

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

#include "FDPS/comm.hpp"
#include "FDPS/comm_for_all_to_all.hpp"
#include "FDPS/time_profile.hpp"
#include "FDPS/util.hpp"
#include "FDPS/ps_macro_defs.h"

namespace ParticleSimulator {

template <S32 DIM>
void SetNumberOfDomainMultiDimension(S32 np[], S32 rank[]) {
  for (S32 i = 0; i < DIMENSION_LIMIT; i++) {
    np[i] = 1;
    rank[i] = 1;
  }
  std::vector<S32> npv;
  npv.resize(DIM);
  S32 np_tmp = Comm::getNumberOfProc();
  for (S32 d = DIM, cid = 0; cid < DIM - 1; d--, cid++) {
    S32 tmp = (S32)pow(np_tmp + 0.000001, (1.0 / d) * 1.000001);
    while (np_tmp % tmp) {
      tmp--;
    }
    npv[cid] = tmp;
    np_tmp /= npv[cid];
  }
  npv[DIM - 1] = np_tmp;
  S32 rank_tmp = Comm::getRank();
  std::sort(npv.begin(), npv.end(), std::greater<S32>());
  for (S32 i = DIM - 1; i >= 0; i--) {
    np[i] = npv[i];
    rank[i] = rank_tmp % np[i];
    rank_tmp /= np[i];
  }
}

class DomainInfo {
 private:
  TimeProfile time_profile_;

  F64vec *pos_sample_tot_;
  F64vec *pos_sample_loc_;

  F64ort *pos_domain_;
  F64ort *pos_domain_temp_;

  S32 *n_smp_array_;
  S32 *n_smp_disp_array_;

  F32 coef_ema_;
  S32 target_number_of_sample_particle_;
  S32 number_of_sample_particle_tot_;
  S32 number_of_sample_particle_loc_;
  S32 n_domain_[DIMENSION_LIMIT];  // in 2-dim, n_domain_[2] is always 1.

  F64ort pos_root_domain_;

  bool first_call_by_initialize;
  bool first_call_by_decomposeDomain;

  S32 boundary_condition_;
  bool periodic_axis_[DIMENSION_LIMIT];  // in 2-dim, periodic_axis_[2] is
                                         // always false.

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  // NEW
  MPI_Comm comm_1d_[DIMENSION_LIMIT];
  MPI_Comm comm_sub_[DIMENSION_LIMIT];
  int rank_1d_[DIMENSION_LIMIT];
  int rank_sub_[DIMENSION_LIMIT];
  int n_proc_sub_[DIMENSION_LIMIT];
#endif

  void calculateBoundaryOfDomainX(const S32 np, const F64vec pos_sample[],
                                  const S32 istart, const S32 iend, F64 &xlow,
                                  F64 &xhigh);

  void calculateBoundaryOfDomainY(const S32 np, const F64vec pos_sample[],
                                  const S32 istart, const S32 iend, F64 &xlow,
                                  F64 &xhigh);

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  void calculateBoundaryOfDomainZ(const S32 np, const F64vec pos_sample[],
                                  const S32 istart, const S32 iend, F64 &xlow,
                                  F64 &xhigh);
#endif

 public:
#ifdef TEST_VARIADIC_TEMPLATE
  void DomainInfoDummyFunc() {}
#endif

  TimeProfile getTimeProfile() const { return time_profile_; }

  void clearTimeProfile() { time_profile_.clear(); }

  DomainInfo();
  ~DomainInfo();

  void initialize(const F32 coef_ema = FDPS_DFLT_VAL_COEF_EMA);

  void setNumberOfDomainMultiDimension(const S32 nx, const S32 ny,
                                       const S32 nz = 1);

  void setDomain(const S32 nx, const S32 ny, const S32 nz = 1) {
    setNumberOfDomainMultiDimension(nx, ny, nz);
  }

  template <class Tpsys>
  void collectSampleParticle(Tpsys &psys, const bool clear, const F32 weight) {
    F64 time_offset = GetWtime();
    if (psys.getFirstCallByDomainInfoCollectSampleParticle()) {
      F64vec *temp_loc = new F64vec[target_number_of_sample_particle_];
      for (S32 i = 0; i < number_of_sample_particle_loc_; i++)
        temp_loc[i] = pos_sample_loc_[i];
      target_number_of_sample_particle_ +=
          psys.getTargetNumberOfSampleParticle();
      delete[] pos_sample_tot_;
      delete[] pos_sample_loc_;

      pos_sample_tot_ = new F64vec[target_number_of_sample_particle_];
      pos_sample_loc_ = new F64vec[target_number_of_sample_particle_];
      for (S32 i = 0; i < number_of_sample_particle_loc_; i++)
        pos_sample_loc_[i] = temp_loc[i];
      delete[] temp_loc;
    }
    if (clear) {
      number_of_sample_particle_loc_ = 0;
    }
    S32 number_of_sample_particle = 0;
    psys.getSampleParticle(number_of_sample_particle,
                           &pos_sample_loc_[number_of_sample_particle_loc_],
                           weight);
    number_of_sample_particle_loc_ += number_of_sample_particle;
    time_profile_.collect_sample_particle += GetWtime() - time_offset;
    return;
  }

  template <class Tpsys>
  void collectSampleParticle(Tpsys &psys, const bool clear) {
    const F32 wgh = psys.getNumberOfParticleLocal();
    collectSampleParticle(psys, clear, wgh);
  }

  template <class Tpsys>
  void collectSampleParticle(Tpsys &psys) {
    const F32 wgh = psys.getNumberOfParticleLocal();
    const bool clear = true;
    collectSampleParticle(psys, clear, wgh);
  }

  // new version multi-dimensional gathering
  void decomposeDomainMultiStep();

  void decomposeDomain();

  template <class Tpsys>
  void decomposeDomainAll(Tpsys &psys, const F32 wgh) {
    const bool clear = true;
    collectSampleParticle(psys, clear, wgh);
    decomposeDomain();
  }

  template <class Tpsys>
  void decomposeDomainAll(Tpsys &psys) {
    const F32 wgh = psys.getNumberOfParticleLocal();
    const bool clear = true;
    collectSampleParticle(psys, clear, wgh);
    decomposeDomain();
  }

  void getRootDomain(FILE *fp) {
    fprintf(fp, "%+e %+e %+e\n", pos_root_domain_.low_[0],
            pos_root_domain_.low_[1], pos_root_domain_.low_[2]);
    fprintf(fp, "%+e %+e %+e\n", pos_root_domain_.high_[0],
            pos_root_domain_.high_[1], pos_root_domain_.high_[2]);

    return;
  }

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
  void getSampleParticleLocal(FILE *fp) {
    for (S32 i = 0; i < number_of_sample_particle_loc_; i++) {
      fprintf(fp, "%+e %+e %+e\n", pos_sample_loc_[i].x, pos_sample_loc_[i].y,
              pos_sample_loc_[i].z);
    }
    return;
  }

  void getSampleParticleTotal(FILE *fp) {
    for (S32 i = 0; i < number_of_sample_particle_tot_; i++) {
      fprintf(fp, "%+e %+e %+e\n", pos_sample_tot_[i].x, pos_sample_tot_[i].y,
              pos_sample_tot_[i].z);
    }
    return;
  }
#endif

  void getPosDomainTotal(FILE *fp) {
    S32 nproc = Comm::getNumberOfProc();
    for (S32 i = 0; i < nproc; i++) {
      for (S32 k = 0; k < DIMENSION; k++)
        fprintf(fp, "%+e ", pos_domain_[i].low_[k]);
      for (S32 k = 0; k < DIMENSION; k++)
        fprintf(fp, "%+e ", pos_domain_[i].high_[k]);
      fprintf(fp, "\n");
    }
  }

  S32 *getPointerOfNDomain() { return n_domain_; };

  F64ort *getPointerOfPosDomain() { return pos_domain_; };

  // A. Tanikawa need this method for Particle Mesh...
  S32 getNDomain(const S32 dim) const { return n_domain_[dim]; };

  // for DEBUG
  F64vec &getPosSample(const S32 id = 0) { return pos_sample_tot_[id]; }
  // for DEBUG // A. Tanikawa would not like to delete this method...
  F64ort &getPosDomain(const S32 id = 0) const { return pos_domain_[id]; }
  // for DEBUG
  void setPosDomain(const S32 id, const F64ort &pos) { pos_domain_[id] = pos; }

  void setBoundaryCondition(enum BOUNDARY_CONDITION bc);

  S32 getBoundaryCondition() const { return boundary_condition_; }

  void setPosRootDomain(const F64vec &low, const F64vec &high);

  const F64ort getPosRootDomain() const { return pos_root_domain_; }

  void getPeriodicAxis(bool pa[]) const {
    for (S32 i = 0; i < DIMENSION; i++) pa[i] = periodic_axis_[i];
  }

  template <class Tpsys>
  bool checkCollectSampleParticleSubset(Tpsys &psys);
  template <class Tpsys>
  bool checkCollectSampleParticleAverage(Tpsys &psys);
  template <class Tpsys>
  bool checkDecomposeDomain(Tpsys &psys);
};

}  // namespace ParticleSimulator
