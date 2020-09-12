#ifndef FDPS_TIME_PROFILE_HPP
#define FDPS_TIME_PROFILE_HPP

#include <iostream>

#include "FDPS/ps_defs.hpp"

namespace ParticleSimulator {

class TimeProfile {
 public:
  F64 collect_sample_particle;
  F64 decompose_domain;
  F64 exchange_particle;
  F64 set_particle_local_tree;
  F64 set_particle_global_tree;
  F64 make_local_tree;
  F64 make_global_tree;
  F64 set_root_cell;
  F64 calc_force;
  F64 calc_moment_local_tree;
  F64 calc_moment_global_tree;
  F64 make_LET_1st;
  F64 make_LET_2nd;
  F64 exchange_LET_1st;
  F64 exchange_LET_2nd;
  F64 write_back;  // new but default

  // not public
  F64 morton_sort_local_tree;
  F64 link_cell_local_tree;
  F64 morton_sort_global_tree;
  F64 link_cell_global_tree;

  F64 make_local_tree_tot;  // make_local_tree + calc_moment_local_tree
  F64 make_global_tree_tot;
  F64 exchange_LET_tot;  // make_LET_1st + make_LET_2nd + exchange_LET_1st +
                         // exchange_LET_2nd

  F64 calc_force__core__walk_tree;
  F64 calc_force__core__keep_list;
  F64 calc_force__core__copy_ep;
  F64 calc_force__core__dispatch;
  F64 calc_force__core__retrieve;

  F64 calc_force__make_ipgroup;
  F64 calc_force__core;
  F64 calc_force__copy_original_order;

  F64 exchange_particle__find_particle;
  F64 exchange_particle__exchange_particle;

  F64 decompose_domain__sort_particle_1st;
  F64 decompose_domain__sort_particle_2nd;
  F64 decompose_domain__sort_particle_3rd;
  F64 decompose_domain__gather_particle;

  F64 decompose_domain__setup;
  F64 decompose_domain__determine_coord_1st;
  F64 decompose_domain__migrae_particle_1st;
  F64 decompose_domain__determine_coord_2nd;
  F64 decompose_domain__determine_coord_3rd;
  F64 decompose_domain__exchange_pos_domain;

  F64 exchange_LET_1st__a2a_n;
  F64 exchange_LET_1st__icomm_sp;
  F64 exchange_LET_1st__a2a_sp;
  F64 exchange_LET_1st__icomm_ep;
  F64 exchange_LET_1st__a2a_ep;

  F64 add_moment_as_sp_local;
  F64 add_moment_as_sp_global;

  TimeProfile();
  TimeProfile operator+(const TimeProfile &rhs) const;
  F64 getTotalTime() const;
  void clear();
  void dump(std::ostream &fout = std::cout) const;
};

}  // namespace ParticleSimulator

#endif  // FDPS_TIME_PROFILE_HPP