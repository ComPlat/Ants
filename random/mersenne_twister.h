/* Portions derived from R: Copyright (C) 2000–2023 The R Foundation */
#ifndef ANTS_RANDOM_H
#define ANTS_RANDOM_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define N 624
#define i2_32m1 2.328306437080797e-10

#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

// Tempering parameters
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)


typedef struct {
  unsigned int i_seed[N + 1];
  int mti;
  bool initial;
} MersenneTwister;

double uniform(int seed, bool init);

#endif // ! ANTS_RANDOM_H

