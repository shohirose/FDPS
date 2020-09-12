#include "FDPS/functions.hpp"

#include "FDPS/comm_table.hpp"
#include "FDPS/memory_pool.hpp"

namespace ParticleSimulator {

void Initialize(int &argc, char **&argv, int64_t memoryPoolSize) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI_Init(&argc, &argv);
#endif
  MemoryPool::initialize(memoryPoolSize);

  if (Comm::getRank() == 0) {
    std::cerr
        << "     //==================================\\\\\n"
        << "     ||                                  ||\n"
        << "     || ::::::: ::::::. ::::::. .::::::. ||\n"
        << "     || ::      ::    : ::    : ::       ||\n"
        << "     || ::::::  ::    : ::::::'  `:::::. ||\n"
        << "     || ::      ::::::' ::      `......' ||\n"
        << "     ||     Framework for Developing     ||\n"
        << "     ||        Particle Simulator        ||\n"
        << "     ||     Version 5.0g (2019/09)       ||\n"
        << "     \\\\==================================//\n\n"
        << "       Home   : https://github.com/fdps/fdps\n"
        << "       E-mail : fdps-support@mail.jmlab.jp\n"
        << "       Licence: MIT (see, "
           "https://github.com/FDPS/FDPS/blob/master/LICENSE)\n"
        << "       Note   : Please cite the following papers.\n"
        << "                - Iwasawa et al. (2016, Publications of the "
           "Astronomical Society of Japan, 68, 54)\n"
        << "                - Namekata et al. (2018, Publications of the "
           "Astronomical Society of Japan, 70, 70)\n\n"
        << "       Copyright (C) 2015\n"
        << "         Masaki Iwasawa, Ataru Tanikawa, Natsuki Hosono,\n"
        << "         Keigo Nitadori, Takayuki Muranushi, Daisuke Namekata,\n"
        << "         Kentaro Nomura, Junichiro Makino and many others\n"
        << "******** FDPS has successfully begun. ********\n"
        << std::endl;
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

}  // namespace ParticleSimulator
