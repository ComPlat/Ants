#ifndef ANTS_ZMAT2XYZ_H
#define ANTS_ZMAT2XYZ_H

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

bool gmetry(int natoms, double* geo_data,
            int* na, int* nb, int* nc,
            double* xyz);

typedef struct {
  int natoms;
  char** symbols;
  int* na;
  int* nb;
  int* nc;
  double* bond_len;
  double* bond_ang;
  double* dihed_ang;
  int* atomic_numbers;
} ZMatrix;

ZMatrix* read_zmatrix(const char* filename);

void free_zmatrix(ZMatrix* zm);

#endif // ! ANTS_ZMAT2XYZ_H
