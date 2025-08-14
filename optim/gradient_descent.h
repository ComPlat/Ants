#ifndef ANTS_GRADIENT_DECENT_H
#define ANTS_GRADIENT_DECENT_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

typedef struct {
  double* parameters;
  double global_best_error;
} ResultLocal;

typedef void (*loss_and_grad_f) (double* parameters, double* error, double* gradient,
                                 bool calc_grad, void* user_data);


int gradient_descent(loss_and_grad_f f, double* x0, int n,
                     double lr, double beta, double tol, int max_iter,
                     void* user_data, Result* res);

#endif // ! ANTS_GRADIENT_DECENT_H
