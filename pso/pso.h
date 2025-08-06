#ifndef ANTS_PSO_H
#define ANTS_PSO_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "mersenne_twister.h"

typedef double (*loss_fct) (double* parameters, void* user_data);

typedef struct {
  double* parameters;
  double global_best_error;
} Result;

void pso(double* lb, double* ub,
         int ngen, int npop, int npar,
         double error_threshold,
         void* user_data, loss_fct lf,
         int seed, Result* res);

#endif // ! ANTS_PSO_H
