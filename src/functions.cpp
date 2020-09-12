#include "FDPS/functions.hpp"

#include <sstream>
#include <string>

#include "FDPS/comm.hpp"
#include "FDPS/memory_pool.hpp"

namespace ParticleSimulator {

void Initialize(int &argc, char **&argv, int64_t memoryPoolSize) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Init(&argc, &argv);
#endif
  MemoryPool::initialize(memoryPoolSize);

  if (Comm::getRank() == 0) {
    PrintLicense(std::cerr);
    std::cerr << "******** FDPS has successfully begun. ********" << std::endl;
  }
}

void Finalize() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Finalize();
#endif
  bool flag_monar = false;
  if (Comm::getRank() == 0) {
    std::cerr << "******** FDPS has successfully finished. ********"
              << std::endl;
  }
}

void PrintLicense(std::ostream &os) {
  std::stringstream ss;
  // clang-format off
  ss << "//==================================\\\\\n"
     << "||                                  ||\n"
     << "|| ::::::: ::::::. ::::::. .::::::. ||\n"
     << "|| ::      ::    : ::    : ::       ||\n"
     << "|| ::::::  ::    : ::::::'  `:::::. ||\n"
     << "|| ::      ::::::' ::      `......' ||\n"
     << "||     Framework for Developing     ||\n"
     << "||        Particle Simulator        ||\n"
     << "||     Version 5.0g (2019/09)       ||\n"
     << "\\\\==================================//\n\n"
     << "Home   : https://github.com/fdps/fdps\n"
     << "E-mail : fdps-support@mail.jmlab.jp\n"
     << "Licence: MIT (https://github.com/FDPS/FDPS/blob/master/LICENSE)\n"
     << "Note   : Please cite the following papers.\n"
     << "  - Iwasawa et al. (2016, Publ. Astron. Soc. Japan, 68, 54)\n"
     << "  - Namekata et al. (2018, Publ. Astron. Soc. Japan, 70, 70)\n\n"
     << "Copyright (C) 2015\n"
     << "  Masaki Iwasawa, Ataru Tanikawa, Natsuki Hosono,\n"
     << "  Keigo Nitadori, Takayuki Muranushi, Daisuke Namekata,\n"
     << "  Kentaro Nomura, Junichiro Makino and many others\n";
  // clang-format on
  os << ss.str();
}

}  // namespace ParticleSimulator
