#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "xtb.h"

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
                            double accuracy, int max_iter, bool verbose);

void destroy_settings_loss_fct(settings_loss_fct* slf);

double loss_fct(settings_loss_fct* slf, const double* coord);

typedef struct {
  double* coord;
  double global_best_error;
} Result;


void pso(double* lb, double* ub,
           int ngen, int npop, int npar,
           double error_threshold,
           const int* atomic_numbers, settings_loss_fct* sf,
           int seed);
