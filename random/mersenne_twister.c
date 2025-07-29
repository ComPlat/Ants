/*
 * This Mersenne Twister RNG implementation is adapted from the R source code.
 * All references to R's C API (e.g., SEXP, PROTECT/UNPROTECT) and integration
 * with R's RNG system have been removed. The remaining code has been refactored
 * for standalone use in pure C, with a simplified interface for initializing and
 * generating uniform random numbers in [0, 1).
 *
 * Original source: R Core Team (https://cran.r-project.org/)
 * License: GPL-2 | GPL-3
*/

#include "mersenne_twister.h"

double fixup(double obj) {
  if (obj == 0.0) {
    return 0.5 * i2_32m1;
  }
  if ((1.0 - obj) <= 0.0) {
    return 1.0 - 0.5 * i2_32m1;
  }
  return obj;
}

int init_scrambling(int seed) {
  for (int i = 0; i < 50; i++) {
    seed = seed * 69069 + 1;
  }
  return seed;
}


void rng_init(MersenneTwister *mt, int seed) {
  for (size_t i = 0; i <= N; i++) {
    seed = seed * 69069 + 1;
    mt->i_seed[i] = seed;
  }
}

void fixup_seeds(MersenneTwister *mt) {
  if (mt->initial) {
    mt->i_seed[0] = 624;
    mt->initial = false;
  }
  if (mt->i_seed[0] <= 0) {
    mt->i_seed[0] = 624;
  }
  bool all_zero = true;
  for (size_t i = 1; i <= N; i++) {
    if (mt->i_seed[i] != 0) {
      all_zero = false;
      break;
    }
  }
  if (all_zero) {
    int new_seed = (int)time(NULL);
    rng_init(mt, new_seed);
  }
}

MersenneTwister mersenne_twister(int seed) {
  MersenneTwister mt;
  mt.initial = true;
  int scrambled_seed = init_scrambling(seed);
  rng_init(&mt, scrambled_seed);
  fixup_seeds(&mt);
  return mt;
}

// Initializing the array with a seed
void MT_sgenrand(unsigned int *mt, int seed, int *mti) {
  int i;
  for (i = 0; i < N; i++) {
    mt[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    mt[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  *mti = N;
}

// Initialization by "sgenrand()" is an example.Theoretically,
//     there are 2 ^
//         19937 - 1 possible states as an initial state.Essential bits in
//                 "seed_array[]" is following 19937 bits
//     : (seed_array[0] & UPPER_MASK),
//     seed_array[1], ...,
//     seed_array[N - 1].(seed_array[0] & LOWER_MASK) is
//     discarded.Theoretically, (seed_array[0] & UPPER_MASK), seed_array[1],
//     ..., seed_array[N - 1] can take any values except all zeros.
double MT_genrand(MersenneTwister *MT) {
  unsigned int y;
  unsigned int *mt = MT->i_seed + 1;
  static unsigned int mag01[2] = {0x0, MATRIX_A};
  int mti = N + 1;
  mti = MT->i_seed[0];

  if (mti >= N) { // generate N words at one time
    int kk;

    if (mti == N + 1) {
      MT_sgenrand(mt, 4357, &mti);
    }
    for (kk = 0; kk < N - M; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (; kk < N - 1; kk++) {
      y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

    mti = 0;
  }

  y = mt[mti++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);
  MT->i_seed[0] = mti;

  return ((double)y * 2.3283064365386963e-10); /* reals: [0,1)-interval */
}

double uniform(int seed, bool init) {
  static MersenneTwister mt;
  if(init) mt = mersenne_twister(seed);
  return MT_genrand(&mt);
}
