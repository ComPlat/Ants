#include "gradient_descent.h"

static bool converged(const double* x_new, const double* x, int n, double tol) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double d = x_new[i] - x[i];
    sum += d * d;
  }
  return sqrt(sum) < tol;
}

static void backtrack(loss_and_grad_f f,
                      const double* x, double* x_new, double* v, int n,
                      double xerror, double* lr, void* user_data) {
  const double beta_bt = 0.5;
  double error = 0.0;

  while (true) {
    for (int i = 0; i < n; i++) x_new[i] = x[i] + v[i];
    f(x_new, &error, NULL, false, user_data);
    if (error < xerror) break;

    *lr *= beta_bt;
    if (*lr < 1e-12) break;
    for (int i = 0; i < n; i++) v[i] *= beta_bt;
  }
}

int gradient_descent(loss_and_grad_f f, double* x0, int n,
                     double lr, double beta, double tol, int max_iter,
                     void* user_data, Result* res) {

  const int min_iter = 3;
  double* x = calloc(n, sizeof(double));
  if (!x) { fprintf(stderr, "alloc x failed\n"); return -1; }
  for (int i = 0; i < n; i++) x[i] = x0[i];

  double* v = calloc(n, sizeof(double));
  if (!v) { free(x); fprintf(stderr, "alloc v failed\n"); return -1; }

  double* grad = calloc(n, sizeof(double));
  if (!grad) { free(x); free(v); fprintf(stderr, "alloc grad failed\n"); return -1; }

  double* x_new = calloc(n, sizeof(double));
  if (!x_new) { free(x); free(v); free(grad); fprintf(stderr, "alloc x_new failed\n"); return -1; }

  double error = 0.0;

  for (int it = 0; it < max_iter; it++) {
    f(x, &error, grad, true, user_data);

    for (int j = 0; j < n; j++)
      v[j] = v[j] * beta - lr * grad[j];

    backtrack(f, x, x_new, v, n, error, &lr, user_data);

    if (converged(x_new, x, n, tol) && it > min_iter) {
      res->parameters = (double*) calloc(n, sizeof(double));
      for (int j = 0; j < n; j++) res->parameters[j] = x_new[j];
      free(x); free(v); free(grad); free(x_new);
      printf("Converged in %d iterations\n", it);
      return 0;
    }

    for (int j = 0; j < n; j++) x[j] = x_new[j];
    if (it % 10 == 0) {
      printf("Iteration: %d, error: %8.6f\n", it, error);
    }
  }

  for (int j = 0; j < n; j++) x0[j] = x[j];
  free(x); free(v); free(grad); free(x_new);
  printf("No convergence\n");
  return 1;
}
