#ifndef ANTS_CALC_ENERGY_H
#define ANTS_CALC_ENERGY_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "xtb.h"
#include "convert.h"

typedef double (*loss_fct) (double* parameters, void* user_data);

typedef struct {
  int natoms;
  int* attyp;
  xtb_TEnvironment env;
  xtb_TMolecule mol;
  xtb_TCalculator calc;
  xtb_TResults res;
  double electronic_temperature;
  double accuracy;
  int max_iter;
  bool verbose;
  int* na;
  int* nb;
  int* nc;
  double* coord;
} settings_loss_fct;

void init_settings_loss_fct(settings_loss_fct* slf, const int natoms, const int* attyp,
                            int* na, int* nb, int* nc,
                            double electronic_temperature,
                            double accuracy, int max_iter, bool verbose);

void destroy_settings_loss_fct(settings_loss_fct* slf);

double calc_energy(double* coord, void* data);

void calc_energy_and_grad(double* coord, double* energy, double* gradient, bool calc_grad, void* data);

#endif // ! ANTS_CALC_ENERGY_H
