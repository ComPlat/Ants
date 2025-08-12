#include "utils.h"

void stats(const double *v, int n, double *m, double *sd) {
  if (n <= 0) { *m = 0.0; *sd = 0.0; return; }

  double best = v[0];
  for (int i = 1; i < n; ++i) if (v[i] < best) best = v[i];
  const double denom = fabs(best) + 1e-12;

  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += (v[i] - best) / denom;
  double mu = sum / (double)n;

  double s2 = 0.0;
  if (n > 1) {
    for (int i = 0; i < n; ++i) {
      double r = (v[i] - best) / denom;
      double d = r - mu;
      s2 += d*d;
    }
    s2 /= (double)(n - 1);
  }
  *m = mu;
  *sd = (n > 1) ? sqrt(s2) : 0.0;
}
