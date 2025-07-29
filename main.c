#include "pso.h"

int main(void) {
  const int natoms = 3;
  const int attyp[3] = {8, 1, 1}; // O, H, H
  const double coord[3 * 3] = {
    0.0000,  0.0000,  0.0000,  // Oxygen
    0.7586,  0.0000,  0.5043,  // Hydrogen 1
   -0.7586,  0.0000,  0.5043   // Hydrogen 2
  };
  double accuracy = 1.0;
  double electronic_temperature = 300.0;
  const int max_iter = 50;
  bool verbose = false;
  settings_loss_fct slf;
  init_settings_loss_fct(&slf, natoms, attyp, electronic_temperature, accuracy, max_iter, verbose);
  double energy = loss_fct(&slf, coord);
  printf("GFN2-xTB single-point energy (H2O): %.10f Eh\n", energy);

  double lb[] = {
    -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0
  };
  double ub[] = {
    1.0, 1.0, 1.0,
    1.0, 1.0, 1.0,
    1.0, 1.0, 1.0
  };

  // instead of x,y,z --> one 3D structure and optimize the length and angle to the adjacent atoms
  // PSO optimizes length and angle of binding; Calculate x,y,z coordinates based on length and angle --> call xtb
  // dont optimize the atoms but the bindings
  pso(lb, ub, 1000, 40, 9, -100, attyp, &slf, 1235);

  destroy_settings_loss_fct(&slf);
  return EXIT_SUCCESS;
}
