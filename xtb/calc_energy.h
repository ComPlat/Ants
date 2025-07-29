#ifndef ANTS_CALC_ENERGY_H
#define ANTS_CALC_ENERGY_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "xtb.h"

typedef double (*loss_fct) (const double* parameters, void* user_data);

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

double calc_energy(const double* coord, void* data);

#endif // ! ANTS_CALC_ENERGY_H
