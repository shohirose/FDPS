#include "FDPS/time_profile.hpp"

namespace ParticleSimulator {

void TimeProfile::dump(std::ostream &fout) const {
  fout << "collect_sample_particle= " << collect_sample_particle << std::endl;
  fout << "decompose_domain= " << decompose_domain << std::endl;
  fout << "exchange_particle= " << exchange_particle << std::endl;
  fout << "set_particle_local_tree= " << set_particle_local_tree << std::endl;
  fout << "set_particle_global_tree= " << set_particle_global_tree << std::endl;
  fout << "make_local_tree= " << make_local_tree << std::endl;
  fout << "make_global_tree= " << make_global_tree << std::endl;
  fout << "set_root_cell= " << set_root_cell << std::endl;
  fout << "calc_force= " << calc_force << std::endl;
  fout << "calc_moment_local_tree= " << calc_moment_local_tree << std::endl;
  fout << "calc_moment_global_tree= " << calc_moment_global_tree << std::endl;
  fout << "make_LET_1st= " << make_LET_1st << std::endl;
  fout << "make_LET_2nd= " << make_LET_2nd << std::endl;
  fout << "exchange_LET_1st= " << exchange_LET_1st << std::endl;
  fout << "exchange_LET_2nd= " << exchange_LET_2nd << std::endl;
  fout << "write_back= " << write_back << std::endl;
}

TimeProfile::TimeProfile() {
  collect_sample_particle = decompose_domain = exchange_particle =
      set_particle_local_tree = set_particle_global_tree = make_local_tree =
          make_global_tree = set_root_cell = calc_force =
              calc_moment_local_tree = calc_moment_global_tree = make_LET_1st =
                  make_LET_2nd = exchange_LET_1st = exchange_LET_2nd = 0.0;
  morton_sort_local_tree = link_cell_local_tree = morton_sort_global_tree =
      link_cell_global_tree = 0.0;
  make_local_tree_tot = make_global_tree_tot = exchange_LET_tot = 0.0;
  calc_force__make_ipgroup = calc_force__core =
      calc_force__copy_original_order = 0.0;

  exchange_particle__find_particle = exchange_particle__exchange_particle = 0.0;

  decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd =
      decompose_domain__sort_particle_3rd = decompose_domain__gather_particle =
          0.0;
  decompose_domain__setup = decompose_domain__determine_coord_1st =
      decompose_domain__migrae_particle_1st =
          decompose_domain__determine_coord_2nd =
              decompose_domain__determine_coord_3rd =
                  decompose_domain__exchange_pos_domain = 0.0;
  exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp =
      exchange_LET_1st__icomm_ep = exchange_LET_1st__icomm_sp =
          exchange_LET_1st__a2a_ep = 0.0;

  add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
  write_back = 0.0;
}

TimeProfile TimeProfile::operator+(const TimeProfile &rhs) const {
  TimeProfile ret;
  ret.collect_sample_particle =
      this->collect_sample_particle + rhs.collect_sample_particle;
  ret.decompose_domain = this->decompose_domain + rhs.decompose_domain;
  ret.exchange_particle = this->exchange_particle + rhs.exchange_particle;
  ret.set_particle_local_tree =
      this->set_particle_local_tree + rhs.set_particle_local_tree;
  ret.set_particle_global_tree =
      this->set_particle_global_tree + rhs.set_particle_global_tree;
  ret.set_root_cell = this->set_root_cell + rhs.set_root_cell;
  ret.make_local_tree = this->make_local_tree + rhs.make_local_tree;
  ret.make_global_tree = this->make_global_tree + rhs.make_global_tree;
  ret.calc_force = this->calc_force + rhs.calc_force;
  ret.calc_moment_local_tree =
      this->calc_moment_local_tree + rhs.calc_moment_local_tree;
  ret.calc_moment_global_tree =
      this->calc_moment_global_tree + rhs.calc_moment_global_tree;
  ret.make_LET_1st = this->make_LET_1st + rhs.make_LET_1st;
  ret.make_LET_2nd = this->make_LET_2nd + rhs.make_LET_2nd;
  ret.exchange_LET_1st = this->exchange_LET_1st + rhs.exchange_LET_1st;
  ret.exchange_LET_2nd = this->exchange_LET_2nd + rhs.exchange_LET_2nd;

  ret.morton_sort_local_tree =
      this->morton_sort_local_tree + rhs.morton_sort_local_tree;
  ret.link_cell_local_tree =
      this->link_cell_local_tree + rhs.link_cell_local_tree;
  ret.morton_sort_global_tree =
      this->morton_sort_global_tree + rhs.morton_sort_global_tree;
  ret.link_cell_global_tree =
      this->link_cell_global_tree + rhs.link_cell_global_tree;

  ret.make_local_tree_tot = this->make_local_tree_tot + rhs.make_local_tree_tot;
  ret.make_global_tree_tot =
      this->make_global_tree_tot + rhs.make_global_tree_tot;
  ret.exchange_LET_tot = this->exchange_LET_tot + rhs.exchange_LET_tot;

  ret.calc_force__core__walk_tree =
      this->calc_force__core__walk_tree + rhs.calc_force__core__walk_tree;
  ret.calc_force__core__keep_list =
      this->calc_force__core__keep_list + rhs.calc_force__core__keep_list;
  ret.calc_force__core__dispatch =
      this->calc_force__core__dispatch + rhs.calc_force__core__dispatch;
  ret.calc_force__core__copy_ep =
      this->calc_force__core__copy_ep + rhs.calc_force__core__copy_ep;
  ret.calc_force__core__retrieve =
      this->calc_force__core__retrieve + rhs.calc_force__core__retrieve;

  ret.calc_force__make_ipgroup =
      this->calc_force__make_ipgroup + rhs.calc_force__make_ipgroup;
  ret.calc_force__core = this->calc_force__core + rhs.calc_force__core;
  ret.calc_force__copy_original_order = this->calc_force__copy_original_order +
                                        rhs.calc_force__copy_original_order;

  ret.exchange_particle__find_particle =
      this->exchange_particle__find_particle +
      rhs.exchange_particle__find_particle;
  ret.exchange_particle__exchange_particle =
      this->exchange_particle__exchange_particle +
      rhs.exchange_particle__exchange_particle;

  ret.decompose_domain__sort_particle_1st =
      this->decompose_domain__sort_particle_1st +
      rhs.decompose_domain__sort_particle_1st;
  ret.decompose_domain__sort_particle_2nd =
      this->decompose_domain__sort_particle_2nd +
      rhs.decompose_domain__sort_particle_2nd;
  ret.decompose_domain__sort_particle_3rd =
      this->decompose_domain__sort_particle_3rd +
      rhs.decompose_domain__sort_particle_3rd;
  ret.decompose_domain__gather_particle =
      this->decompose_domain__gather_particle +
      rhs.decompose_domain__gather_particle;

  ret.decompose_domain__setup =
      this->decompose_domain__setup + rhs.decompose_domain__setup;
  ret.decompose_domain__determine_coord_1st =
      this->decompose_domain__determine_coord_1st +
      rhs.decompose_domain__determine_coord_1st;
  ret.decompose_domain__migrae_particle_1st =
      this->decompose_domain__migrae_particle_1st +
      rhs.decompose_domain__migrae_particle_1st;
  ret.decompose_domain__determine_coord_2nd =
      this->decompose_domain__determine_coord_2nd +
      rhs.decompose_domain__determine_coord_2nd;
  ret.decompose_domain__determine_coord_3rd =
      this->decompose_domain__determine_coord_3rd +
      rhs.decompose_domain__determine_coord_3rd;
  ret.decompose_domain__exchange_pos_domain =
      this->decompose_domain__exchange_pos_domain +
      rhs.decompose_domain__exchange_pos_domain;

  ret.exchange_LET_1st__a2a_n =
      this->exchange_LET_1st__a2a_n + rhs.exchange_LET_1st__a2a_n;
  ret.exchange_LET_1st__icomm_ep =
      this->exchange_LET_1st__icomm_ep + rhs.exchange_LET_1st__icomm_ep;
  ret.exchange_LET_1st__a2a_sp =
      this->exchange_LET_1st__a2a_sp + rhs.exchange_LET_1st__a2a_sp;
  ret.exchange_LET_1st__icomm_sp =
      this->exchange_LET_1st__icomm_sp + rhs.exchange_LET_1st__icomm_sp;
  ret.exchange_LET_1st__a2a_ep =
      this->exchange_LET_1st__a2a_ep + rhs.exchange_LET_1st__a2a_ep;

  ret.add_moment_as_sp_local =
      this->add_moment_as_sp_local + rhs.add_moment_as_sp_local;
  ret.add_moment_as_sp_global =
      this->add_moment_as_sp_global + rhs.add_moment_as_sp_global;

  ret.write_back = this->write_back + rhs.write_back;
  return ret;
}

F64 TimeProfile::getTotalTime() const {
  return collect_sample_particle + decompose_domain + exchange_particle +
         set_particle_local_tree + set_particle_global_tree + make_local_tree +
         make_global_tree + set_root_cell + calc_force +
         calc_moment_local_tree + calc_moment_global_tree + make_LET_1st +
         make_LET_2nd + exchange_LET_1st + exchange_LET_2nd +
         morton_sort_local_tree + link_cell_local_tree +
         morton_sort_global_tree + link_cell_global_tree +
         add_moment_as_sp_local + add_moment_as_sp_global + write_back;
}

void TimeProfile::clear() {
  collect_sample_particle = decompose_domain = exchange_particle =
      make_local_tree = make_global_tree = set_particle_local_tree =
          set_particle_global_tree = set_root_cell = calc_force =
              calc_moment_local_tree = calc_moment_global_tree = make_LET_1st =
                  make_LET_2nd = exchange_LET_1st = exchange_LET_2nd = 0.0;
  morton_sort_local_tree = link_cell_local_tree = morton_sort_global_tree =
      link_cell_global_tree = 0.0;
  make_local_tree_tot = make_global_tree_tot = exchange_LET_tot = 0.0;
  calc_force__core__walk_tree = 0.0;
  calc_force__core__keep_list = 0.0;
  calc_force__core__copy_ep = 0.0;
  calc_force__core__dispatch = 0.0;
  calc_force__core__retrieve = 0.0;
  calc_force__make_ipgroup = calc_force__core =
      calc_force__copy_original_order = 0.0;
  exchange_particle__find_particle = exchange_particle__exchange_particle = 0.0;

  decompose_domain__sort_particle_1st = decompose_domain__sort_particle_2nd =
      decompose_domain__sort_particle_3rd = decompose_domain__gather_particle =
          0.0;
  decompose_domain__setup = decompose_domain__determine_coord_1st =
      decompose_domain__migrae_particle_1st =
          decompose_domain__determine_coord_2nd =
              decompose_domain__determine_coord_3rd =
                  decompose_domain__exchange_pos_domain = 0.0;
  exchange_LET_1st__a2a_n = exchange_LET_1st__a2a_sp =
      exchange_LET_1st__icomm_ep = exchange_LET_1st__icomm_sp =
          exchange_LET_1st__a2a_ep = 0.0;

  add_moment_as_sp_local = add_moment_as_sp_global = 0.0;
  write_back = 0.0;
}

}  // namespace ParticleSimulator
