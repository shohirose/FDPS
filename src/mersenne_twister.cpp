#include "FDPS/MT.hpp"

namespace ParticleSimulator {

void MT::init_by_array(unsigned long init_key[], int key_length) {
  int i, j, k;
  init_genrand(19650218UL);
  i = 1;
  j = 0;
  k = (MT_N > key_length ? MT_N : key_length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) +
            init_key[j] + j;  // non linear
    mt[i] &= 0xffffffffUL;    // for WORDSIZE > 32 machines
    i++;
    j++;
    if (i >= MT_N) {
      mt[0] = mt[MT_N - 1];
      i = 1;
    }
    if (j >= key_length) j = 0;
  }
  for (k = MT_N - 1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) -
            i;              // non linear
    mt[i] &= 0xffffffffUL;  // for WORDSIZE > 32 machines
    i++;
    if (i >= MT_N) {
      mt[0] = mt[MT_N - 1];
      i = 1;
    }
  }
  mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array
}

unsigned long MT::genrand_int32() {
  unsigned long y;
  static unsigned long mag01[2] = {0x0UL, MATRIX_A};
  // mag01[x] = x * MATRIX_A  for x=0,1

  if (mti >= MT_N) {  // generate N words at one time
    int kk;

    if (mti == MT_N + 1)     // if init_genrand() has not been called,
      init_genrand(5489UL);  // a default initial seed is used

    for (kk = 0; kk < MT_N - MT_M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < MT_N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
  }

  y = mt[mti++];

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

void MT::init_genrand(unsigned long s) {
  unsigned long *tmp_mt = getInstance().mt;
  int &tmp_mti = getInstance().mti;

  tmp_mt[0] = s & 0xffffffffUL;
  for (tmp_mti = 1; tmp_mti < MT_N; tmp_mti++) {
    tmp_mt[tmp_mti] =
        (1812433253UL * (tmp_mt[tmp_mti - 1] ^ (tmp_mt[tmp_mti - 1] >> 30)) +
         tmp_mti);
    tmp_mt[tmp_mti] &= 0xffffffffUL;
  }
}

void MTTS::init_by_array(unsigned long init_key[], int key_length) {
  int i, j, k;
  init_genrand(19650218UL);
  i = 1;
  j = 0;
  k = (MT_N > key_length ? MT_N : key_length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) +
            init_key[j] + j;  // non linear
    mt[i] &= 0xffffffffUL;    // for WORDSIZE > 32 machines
    i++;
    j++;
    if (i >= MT_N) {
      mt[0] = mt[MT_N - 1];
      i = 1;
    }
    if (j >= key_length) j = 0;
  }
  for (k = MT_N - 1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) -
            i;              // non linear
    mt[i] &= 0xffffffffUL;  // for WORDSIZE > 32 machines /
    i++;
    if (i >= MT_N) {
      mt[0] = mt[MT_N - 1];
      i = 1;
    }
  }
  mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array
}

unsigned long MTTS::genrand_int32() {
  unsigned long y;
  static unsigned long mag01[2] = {0x0UL, MATRIX_A};
  // mag01[x] = x * MATRIX_A  for x=0,1

  if (mti >= MT_N) {  // generate N words at one time
    int kk;

    if (mti == MT_N + 1)     // if init_genrand() has not been called,
      init_genrand(5489UL);  // a default initial seed is used

    for (kk = 0; kk < MT_N - MT_M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < MT_N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mti = 0;
  }

  y = mt[mti++];

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

void MTTS::init_genrand(unsigned long s) {
  unsigned long *tmp_mt = mt;
  int &tmp_mti = mti;

  tmp_mt[0] = s & 0xffffffffUL;
  for (tmp_mti = 1; tmp_mti < MT_N; tmp_mti++) {
    tmp_mt[tmp_mti] =
        (1812433253UL * (tmp_mt[tmp_mti - 1] ^ (tmp_mt[tmp_mti - 1] >> 30)) +
         tmp_mti);
    tmp_mt[tmp_mti] &= 0xffffffffUL;
  }
}

}  // namespace ParticleSimulator
