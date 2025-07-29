#include "pso.h"
#include "calc_energy.h"

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
  void* data = &slf;
  double energy = calc_energy(coord, data);
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

  Result res = {0};
  pso(lb, ub, 1000, 40, 9, -100, data, calc_energy, 1235, &res);

  printf("[");
  for (int i = 0; i < (natoms*3); i++) {
    printf("%0.1f, ", res.parameters[i]);
  }
  printf("]\n");

  destroy_settings_loss_fct(&slf);
  return EXIT_SUCCESS;
}
