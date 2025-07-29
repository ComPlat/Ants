#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "xtb.h"

static double POS_INF = 1.0 /0.0;

typedef struct {
  int natoms;
  const int* attyp;
  xtb_TEnvironment env;
  xtb_TMolecule mol;
  xtb_TCalculator calc;
  xtb_TResults res;
  double electronic_temperature;
  double accuracy;
  int max_iter;
  bool verbose;
} settings_loss_fct;

void init_settings_loss_fct(settings_loss_fct* slf, const int natoms,
                            const int* attyp, double electronic_temperature,
                            double accuracy, int max_iter, bool verbose) {
  slf->natoms = natoms;
  slf->attyp = attyp;
  slf->env = xtb_newEnvironment();
  slf->calc = xtb_newCalculator();
  slf->res = xtb_newResults();
  slf->electronic_temperature = electronic_temperature;
  slf->accuracy = accuracy;
  slf->max_iter = max_iter;
  slf->verbose = verbose;
  if (!slf->verbose) {
    xtb_setVerbosity(slf->env, XTB_VERBOSITY_MUTED);
  }
  slf->mol = NULL;
}

void destroy_settings_loss_fct(settings_loss_fct* slf) {
  xtb_delete(slf->res);
  xtb_delete(slf->calc);
  xtb_delete(slf->mol);
  xtb_delete(slf->env);
}

double loss_fct(settings_loss_fct* slf, const double* coord) {
  if (slf->mol == NULL) {
    slf->mol = xtb_newMolecule(slf->env, &slf->natoms, slf->attyp, coord, NULL, NULL, NULL, NULL);
  } else {
    xtb_updateMolecule(slf->env, slf->mol, coord, NULL);
  }
  xtb_loadGFN2xTB(slf->env, slf->mol, slf->calc, NULL);
  xtb_setAccuracy(slf->env, slf->calc, slf->accuracy);
  xtb_setElectronicTemp(slf->env, slf->calc, slf->electronic_temperature);
  xtb_setMaxIter(slf->env, slf->calc, slf->max_iter);
  xtb_singlepoint(slf->env, slf->mol, slf->calc, slf->res);
  if(xtb_checkEnvironment(slf->env)) {
    if (slf->verbose) {
      char buffer[512];
      int bufsize = 512;
      xtb_getError(slf->env, buffer, &bufsize);
      fprintf(stderr, "xTB Error: %s\n", buffer);
    }
    return POS_INF;
  }
  double energy;
  xtb_getEnergy(slf->env, slf->res, &energy);
  return energy;
}

typedef struct {
  double* coord;
  double global_best_error;
} Result;

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

void destroy_swarm(Swarm* s);

Swarm* init_empty_swarm() {
  Swarm* s = (Swarm*) malloc(sizeof(Swarm));
  if (!s) return NULL;
  s->swarm         = NULL;
  s->swarm_bests_params = NULL;
  s->v             = NULL;
  s->swarm_bests   = NULL;
  s->swarm_errors  = NULL;
  s->neighborhood  = NULL;
  return s;
}

void check_allocation(Swarm*s) {
  if (!s || !s->swarm || !s->v || !s->swarm_bests_params || !s->swarm_bests || !s->swarm_errors || !s->neighborhood) {
    fprintf(stderr, "Failed to allocate swarm");
    destroy_swarm(s);
  }
}

Swarm* init_swarm(int npop, int npar) {
  Swarm* s = init_empty_swarm();
  s->npop = npop;
  s->npar = npar;
  s->size = npop * npar;
  s->K = 3;
  s -> size_neighberhood = s->K*s->npop;

  s->swarm         = (double*) malloc(s->size * sizeof(double));
  s->swarm_bests_params = (double*) malloc(s->size * sizeof(double));
  s->v             = (double*) calloc(s->size, sizeof(double));
  s->swarm_bests   = (double*) calloc(s->npop, sizeof(double));
  s->swarm_errors  = (double*) calloc(s->npop, sizeof(double));
  s->neighborhood  = (int*)    calloc(s->size_neighberhood, sizeof(int));

  check_allocation(s);

  return s;
}

void init_parameters(double* lb, double* ub, Swarm* s) {
  for (int i = 0; i < s->npop; i++) {
    for (int j = 0; j < s->npar; j++) {
      double r = (double) rand() / RAND_MAX;
      s->swarm[(i * s->npar) + j] = lb[j] + r * (ub[j] - lb[j]);
      s->swarm_bests_params[(i * s->npar) + j] = s->swarm[(i*s->npar) + j];
    }
  }
}

void eval_swarm(Swarm* s, settings_loss_fct* sf) {
  for (int i = 0; i < s->npop; i++) {
    double* coord = &s->swarm[i * s->npar];
    s->swarm_errors[i] = loss_fct(sf, coord);
  }
}

void calc_neighberhood(Swarm*s) {
  for (int i = 0; i < s->size_neighberhood; i++) s->neighborhood[i] = -1;
  for (int i = 0; i < s->npop; i++) {
    int n_neighbors = rand() % s->K;
    for (int j = 0; j < n_neighbors; j++) {
      // TODO: duplicates are wrong
      s->neighborhood[(i * s->K) + j] = rand() % s->npop;
    }
  }
}

void correct_parameters(double* lb, double* ub, Swarm* s) {
  for (int i = 0; i < s->npop; i++) {
    for (int j = 0; j < s->npar; j++) {
      if (s->swarm[(i * s->npar) + j] < lb[j]) s->swarm[(i * s->npar) + j] = lb[j];
      if (s->swarm[(i * s->npar) + j] > ub[j]) s->swarm[(i * s->npar) + j] = ub[j];
    }
  }
}

int find_best_particle(Swarm*s) {
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
int find_best_particle_in_neighborhood(Swarm*s, int current_particle, int* neighborhood) {
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
    if (error_best_particle < s->swarm_errors[neighborhood[indices[i]]]) {
      error_best_particle = s->swarm_errors[neighborhood[indices[i]]];
      best_particle = neighborhood[indices[i]];
    }
  }
  free(indices);
  return best_particle;
}

void update_swarm(Swarm*s, double w, double cog, double soc) {
  for (int i = 0; i < s->npop; i++) {
    int* current_neighborhood = &s->neighborhood[i * s->K];
    int best_particle = find_best_particle_in_neighborhood(s, i, current_neighborhood);
    double* local_best = &s->swarm[best_particle * s->npar];
    double* swarm_bests_params = &s->swarm_bests_params[i * s->npar];
    double* v_i = &s->v[i * s->npar];
    double* swarm_i = &s->swarm[i * s->npar];
    double r1 = ((double) rand() / RAND_MAX)*0.25;
    double r2 = ((double) rand() / RAND_MAX)*0.25;
    for (int j = 0; j <s->npar; j++) {
      v_i[j] = w * v_i[j] + cog * r1 * (swarm_bests_params[j] - swarm_i[j]) + soc * r2 * (local_best[j] - swarm_i[j]);
      swarm_i[j] = swarm_i[j] + v_i[j];
    }
  }
}

void destroy_swarm(Swarm* s) {
  free(s->swarm); s->swarm = NULL;
  free(s->v); s->v = NULL;
  free(s->swarm_bests_params); s->swarm_bests_params = NULL;
  free(s->swarm_bests); s->swarm_bests = NULL;
  free(s->swarm_errors); s->swarm_errors = NULL;
  free(s->neighborhood); s->neighborhood = NULL;
  free(s); s = NULL;
}

void write_xyz(const char* filename, const int* attyp, const double* coords, int natoms) {
    FILE* f = fopen(filename, "w");
    if (!f) return;
    fprintf(f, "%d\n", natoms);
    fprintf(f, "Generated by PSO\n");
    for (int i = 0; i < natoms; i++) {
        int z = attyp[i];
        double x = coords[3 * i + 0];
        double y = coords[3 * i + 1];
        double zc = coords[3 * i + 2];
        fprintf(f, "%d %.10f %.10f %.10f\n", z, x, y, zc);
    }
    fclose(f);
}

void pso(double* lb, double* ub,
         int ngen, int npop, int npar,
         double error_threshold,
         const int* atomic_numbers, settings_loss_fct* sf,
         int seed) {
  srand(seed);
  Swarm* s= init_swarm(npop, npar);
  if (!s) return;
  init_parameters(lb, ub, s);
  eval_swarm(s, sf);
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
  double* global_best_vec= &s->swarm[global_best * s->npar];
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
    update_swarm(s, w, cog, soc);
    correct_parameters(lb, ub, s);
    eval_swarm(s, sf);
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
      for (int j = 0; j < s->npar; j++) {
        global_best_vec[j] = s->swarm[global_best * s->npar + j];
      }
    } else {
      no_improvement++;
    }
    iter++;
    printf("%.10f  %0.10f\n", global_best_error, s->swarm_errors[global_best]);
  }
    printf("Best Parameters: [");
    double* coord = &s->swarm[global_best*s->npar];
    for (int j = 0; j < s->npar; j++) {
      printf("%.5f", coord[j]);
      if (j < s->npar - 1) printf(", ");
    }
    printf("]\n");

  // TODO: user should pass file name
  write_xyz("best_structure.xyz", sf->attyp, coord, sf->natoms);

  destroy_swarm(s);

  // TODO: return Result struct
}
