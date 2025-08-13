#include "pso.h"

static void stats(const double *v, int n, double *m, double *sd) {
  if (n <= 0) { *m = 0.0; *sd = 0.0; return; }

  double best = v[0];
  for (int i = 1; i < n; i++) if (v[i] < best) best = v[i];
  const double denom = fabs(best) + 1e-12; // if best is 0.0

  double sum = 0.0;
  for (int i = 0; i < n; i++) sum += (v[i] - best) / denom;
  const double mean = sum / (double) n;

  double s = 0.0;
  if (n > 1) {
    for (int i = 0; i < n; i++) {
      double r = (v[i] - best) / denom;
      double d = r - mean;
      s += d*d;
    }
    s /= (double)(n - 1);
  }
  *m = mean;
  *sd = (n > 1) ? sqrt(s) : 0.0;
}

typedef struct {
  int npop;
  int npar;
  int size;
  int K;
  int size_neighberhood;
  double* swarm;
  double* v;
  double* swarm_bests_params;
  double* swarm_bests;
  double* swarm_errors;
  int* neighborhood;
} Swarm;

static void destroy_swarm(Swarm* s);

static Swarm* init_empty_swarm() {
  Swarm* s = (Swarm*) malloc(sizeof(Swarm));
  s->swarm         = NULL;
  s->swarm_bests_params = NULL;
  s->v             = NULL;
  s->swarm_bests   = NULL;
  s->swarm_errors  = NULL;
  s->neighborhood  = NULL;
  return s;
}

static bool check_allocation(Swarm*s) {
  if (!s || !s->swarm || !s->v || !s->swarm_bests_params ||
    !s->swarm_bests || !s->swarm_errors || !s->neighborhood) {
    fprintf(stderr, "Failed to allocate swarm");
    return false;
  }
  return true;
}

static Swarm* init_swarm(int npop, int npar) {
  Swarm* s = init_empty_swarm();
  s->npop = npop;
  s->npar = npar;
  s->size = npop * npar;
  s->K = 4;
  s -> size_neighberhood = s->K*s->npop;

  s->swarm         = (double*) malloc(s->size * sizeof(double));
  s->swarm_bests_params = (double*) malloc(s->size * sizeof(double));
  s->v             = (double*) calloc(s->size, sizeof(double));
  s->swarm_bests   = (double*) calloc(s->npop, sizeof(double));
  s->swarm_errors  = (double*) calloc(s->npop, sizeof(double));
  s->neighborhood  = (int*)    calloc(s->size_neighberhood, sizeof(int));

  if (!check_allocation(s)) {
    destroy_swarm(s);
    return NULL;
  }

  return s;
}

void init_parameters(double* lb, double* ub, Swarm* s) {
  for (int i = 0; i < s->npop; i++) {
    for (int j = 0; j < s->npar; j++) {
      double r = uniform(-1, false);
      s->swarm[(i * s->npar) + j] = lb[j] + r * (ub[j] - lb[j]);
      s->swarm_bests_params[(i * s->npar) + j] = s->swarm[(i*s->npar) + j];
    }
  }
}

static void eval_swarm(Swarm* s, void* user_data, loss_fct lf) {
  for (int i = 0; i < s->npop; i++) {
    double* parameter = &s->swarm[i * s->npar];
    s->swarm_errors[i] = lf(parameter, user_data);
  }
}

static void calc_neighberhood(Swarm*s) {
  for (int i = 0; i < s->size_neighberhood; i++) s->neighborhood[i] = -1;
  for (int i = 0; i < s->npop; i++) {
    int n_neighbors = rand_int(s->K, uniform);
    int* out = &s->neighborhood[i * s->K];
    sample_ints(s->npop, n_neighbors, false, uniform, out);
  }
}

static int find_best_particle(Swarm*s) {
  double error_best_particle = s->swarm_errors[0];
  int best_particle = 0;
  for (int i = 1; i < s->npop; i++) {
    if (s->swarm_errors[i] < error_best_particle) {
      error_best_particle = s->swarm_errors[i];
      best_particle = i;
    }
  }
  return best_particle;
}
static int find_best_particle_in_neighborhood(Swarm*s, int current_particle, int* neighborhood) {
  int counter = 0;
  int* indices = (int*) malloc(s->K * sizeof(int));
  for (int i = 0; i < s->K; i++) indices[i] = -1;
  for (int i = 0; i < s->K; i++) {
    if (neighborhood[i] != -1) {
      counter++;
      indices[i] = i;
      continue;
    }
  }
  if (counter == 0) {
    free(indices);
    return current_particle;
  }

  double error_best_particle = s->swarm_errors[neighborhood[indices[0]]];
  int best_particle = neighborhood[indices[0]];
  for (int i = 1; i < counter; i++) {
    if (error_best_particle > s->swarm_errors[neighborhood[indices[i]]]) {
      error_best_particle = s->swarm_errors[neighborhood[indices[i]]];
      best_particle = neighborhood[indices[i]];
    }
  }
  free(indices);
  return best_particle;
}

static void update_and_eval_swarm(Swarm* s, double w, double cog, double soc,
                       double* lb, double* ub,
                       loss_fct lf, void* user_data) {
  for (int i = 0; i < s->npop; i++) {
    int* current_neighborhood = &s->neighborhood[i * s->K];
    int best_particle = find_best_particle_in_neighborhood(s, i, current_neighborhood);

    double* local_best = &s->swarm_bests_params[best_particle * s->npar];
    double* personal_best = &s->swarm_bests_params[i * s->npar];
    double* v_i = &s->v[i * s->npar];
    double* x_i = &s->swarm[i * s->npar];

    double backup_x[s->npar];
    double backup_v[s->npar];
    for (int j = 0; j < s->npar; j++) {
      backup_x[j] = x_i[j];
      backup_v[j] = v_i[j];
    }

    for (int j = 0; j < s->npar; j++) {
      const double r1 = uniform(-1, false);
      const double r2 = uniform(-1, false);

      v_i[j] = w * v_i[j] + cog * r1 * (personal_best[j] - x_i[j]) + soc * r2 * (local_best[j] - x_i[j]);

      const double k = 0.2;
      const double range = fabs(ub[j] - lb[j]);
      const double pull = fabs(personal_best[j]-x_i[j]) + fabs(local_best[j]-x_i[j]);
      double vmax = k * pull;
      const double vmax_cap = k * range;
      if (vmax > vmax_cap) vmax = vmax_cap;
      double vmax_floor = fmax(1e-12, 1e-6 * range);
      if (vmax < vmax_floor) vmax = vmax_floor;
      if (v_i[j] >  vmax) v_i[j] =  vmax;
      if (v_i[j] < -vmax) v_i[j] = -vmax;
      x_i[j] += v_i[j];

      // Clamp
      if (x_i[j] < lb[j]) { x_i[j] = lb[j] + (lb[j]-x_i[j]); v_i[j] = -0.5*v_i[j]; }
      if (x_i[j] > ub[j]) { x_i[j] = ub[j] - (x_i[j]-ub[j]); v_i[j] = -0.5*v_i[j]; }
      if (x_i[j] < lb[j]) x_i[j] = lb[j];
      if (x_i[j] > ub[j]) x_i[j] = ub[j];
    }

    double err = lf(x_i, user_data);
    if (isinf(err) || isnan(err)) {
      // Restore previous state
      for (int j = 0; j < s->npar; j++) {
        x_i[j] = backup_x[j];
        v_i[j] = backup_v[j];
      }
    } else {
      s->swarm_errors[i] = err;
    }
  }
}

typedef struct {
    int idx;
    double fit;
} RankItem;

static int compare_rankitem(const void *a, const void *b) {
    const RankItem *ra = (const RankItem*)a;
    const RankItem *rb = (const RankItem*)b;
    // Handle NaNs: treat NaN as worst
    const int a_nan = isnan(ra->fit), b_nan = isnan(rb->fit);
    if (a_nan && !b_nan) return 1;
    if (!a_nan && b_nan) return -1;
    // ascending (lower = better)
    if (ra->fit < rb->fit) return -1;
    if (ra->fit > rb->fit) return 1;
    return 0;
}

void refresh_swarm(Swarm* s,
                   void* user_data, loss_fct lf,
                   const double* lb, const double* ub) {
  RankItem *tmp = (RankItem*) malloc(s->npop * sizeof(RankItem));
  for (int i = 0; i < s->npop; i++) {
    tmp[i].idx = i; tmp[i].fit = s->swarm_bests[i];
  }
  qsort(tmp, s->npop, sizeof(RankItem), compare_rankitem);
  int start = (int)floor(s->npop * 0.5);
  for (int t = start; t < s->npop; ++t) {
    int i = tmp[t].idx;
    // reinit particle i
    for (int j = 0; j < s->npar; j++) {
      double r = uniform(-1, false);
      s->swarm[i*s->npar + j] = lb[j] + r * (ub[j] - lb[j]);
      s->v[i*s->npar + j] = 0.01 * (ub[j]-lb[j]) * (uniform(-1,false)-0.5);
      s->swarm_bests_params[i*s->npar + j] = s->swarm[i*s->npar + j];
    }
    // eval & set bests
    double* p = &s->swarm[i * s->npar];
    s->swarm_errors[i] = lf(p, user_data);
    s->swarm_bests[i]  = s->swarm_errors[i];
  }
  free(tmp);
}

static void destroy_swarm(Swarm* s) {
  if (!s) return;
  free(s->swarm); s->swarm = NULL;
  free(s->v); s->v = NULL;
  free(s->swarm_bests_params); s->swarm_bests_params = NULL;
  free(s->swarm_bests); s->swarm_bests = NULL;
  free(s->swarm_errors); s->swarm_errors = NULL;
  free(s->neighborhood); s->neighborhood = NULL;
  free(s); s = NULL;
}

static void fill_result(Result* res, double* best_parameters,
                 double global_best_error, int npar) {
  if(!res->parameters) res->parameters = (double*) malloc(sizeof(double)*npar);
  for (int i = 0; i < npar; i++) res->parameters[i] = best_parameters[i];
  res->global_best_error = global_best_error;
}

void pso(double* lb, double* ub,
         int ngen, int npop, int npar,
         double error_threshold,
         void* user_data, loss_fct lf,
         int seed, Result* res) {
  uniform(seed, true);
  Swarm* s = init_swarm(npop, npar);
  if (!s) {
    return;
  }

  init_parameters(lb, ub, s);
  eval_swarm(s, user_data, lf);
  for (int i = 0; i < s->npop; i++) s->swarm_bests[i] = s->swarm_errors[i];

  const double initial_cog = 2.5;
  const double final_cog = 0.5;
  const double initial_soc = 0.5;
  const double final_soc = 2.5;
  double w = 0.5;
  const double w_max = 0.9;
  const double w_min = 0.4;
  double cog = 0.0;
  double soc = 0.0;

  int global_best = find_best_particle(s);
  double global_best_error = s->swarm_bests[global_best];

  int no_improvement = 0;
  int iter = 1;

  while (iter < ngen) {
    if ((iter == 1) || (no_improvement != 0)) {
      calc_neighberhood(s);
    }
    w = w_max - iter * (w_max - w_min) / ngen;
    cog = initial_cog - (initial_cog - final_cog) * (iter + 1) / ngen;
    soc = initial_soc - (initial_soc - final_soc) * (iter + 1) / ngen;
    update_and_eval_swarm(s, w, cog, soc, lb, ub, lf, user_data);
    for (int i = 0; i < s->npop; i++) {
      if (s->swarm_errors[i] < s->swarm_bests[i]) {
        s->swarm_bests[i] = s->swarm_errors[i];
        for (int j = 0; j < s->npar; j++) {
          s->swarm_bests_params[i * s->npar + j] = s->swarm[i * s->npar + j];
        }
      }
    }
    int new_global_best = find_best_particle(s);
    if (s->swarm_errors[new_global_best] < global_best_error) {
      global_best_error = s->swarm_errors[new_global_best];
      global_best = new_global_best;
      no_improvement = 0;
    } else {
      no_improvement++;
    }
    iter++;
    double mean, sd;
    stats(s->swarm_bests, npop, &mean, &sd);
    printf("Generation: %d/%d, global best: %.3f; rel mean/sd: %.3f%% / %.3f%%; Ratio: %.3f \n",
           iter, ngen, global_best_error, 100.0*mean, 100.0*sd, (mean / sd));
    if (global_best_error <= error_threshold) break;
    double ratio = fabs(mean) / sd;
    if ((ratio > 2.5) || (no_improvement == 10)) {
      printf("Refresh part of the swarm \n");
      refresh_swarm(s, user_data, lf, lb, ub);
    }
  }

  double* best_params = &s->swarm[global_best*s->npar];
  fill_result(res, best_params, global_best_error, npar);
  destroy_swarm(s);
}
