#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "pso.h"

double sphere(double* x, void* data) {
  int n = *(int*)data;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += x[i] * x[i];
  return sum;
}

int main() {
  int npar = 5;
  double lb[5], ub[5];
  for (int i = 0; i < npar; ++i) {
    lb[i] = -5.0;
    ub[i] = 5.0;
  }

  Result result = {0};
  void* data = &npar;
  pso(lb, ub, 1000, 40, npar, 1e-4, data, sphere, 1234, &result);

  double final_error = sphere(result.parameters, &npar);
  printf("Final error: %.8f\n", final_error);
  assert(final_error < 1e-3);
  for (int i = 0; i < npar; ++i) {
    printf("x[%d] = %.6f\n", i, result.parameters[i]);
  }
  free(result.parameters);
  return 0;
}
