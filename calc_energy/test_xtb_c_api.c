#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "xtb.h"
static double POS_INF = 1.0 / 0.0;

typedef struct {
  int natoms;
  const int* attyp;
  const double* coord;
  double charge;
  int uhf;
} Data;

Data init_data(int natoms, const int* attyp,
               const double* coord, double charge, int uhf) {
  Data data = {0};
  data.natoms = natoms;
  data.attyp = attyp;
  data.coord = coord;
  data.charge = charge;
  data.uhf = uhf;
  return data;
}

double calc_energy_and_gradient(Data* data, double** out_gradient) {
  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  double energy = POS_INF;

  xtb_TMolecule mol = xtb_newMolecule(
    env, &(data->natoms), data->attyp, data->coord, &(data->charge), &(data->uhf), NULL, NULL);

  if (xtb_checkEnvironment(env)) goto cleanup;

  xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
  if (xtb_checkEnvironment(env)) goto cleanup;

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env)) goto cleanup;

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env)) goto cleanup;

  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env)) goto cleanup;

  *out_gradient = malloc(3 * data->natoms * sizeof(double));
  xtb_getGradient(env, res, *out_gradient);
  if (xtb_checkEnvironment(env)) goto cleanup;

cleanup:
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    energy = POS_INF;
    if (*out_gradient != NULL) {
      free(*out_gradient);
      *out_gradient = NULL;
    }
  }

  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);

  return energy;
}

int main(int argc, char** argv) {
  int    const natoms = 18;
  int    const attyp[18] = {6, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1, 1};
  double const charge = 0.0;
  int    const uhf = 0;
  double const coord[3*18] = {
    0.0000000000, 0.0000000000, 0.0000000000,
    2.9322094206, 0.0000000000, 0.0000000000,
   -1.1742394969, 2.6864456722, 0.0000000000,
   -0.6051014872, -1.0919282327, 1.6924507605,
   -0.7595265229, -0.9882521577, -1.6896416115,
    0.3479268993, 4.4070164168, 1.7406633633,
   -3.1856271240, 2.5656796203, 0.5541491550,
   -1.0942653918, 3.5070543555, -1.9341932218,
    3.1071060664, 4.6862491743, 0.7881920978,
    0.3334160410, 3.4840279798, 3.6273863886,
   -0.5169672863, 6.3093502430, 1.9575489272,
    3.9651154386, 2.5172269236, -0.9824517572,
    4.3284345078, 4.7931151837, 2.4931262597,
    3.3636550613, 6.4705907998, -0.2880493804,
    3.1794030178, 2.8029736653, -2.9133563500,
    6.0567576610, 2.5147967551, -1.1393672182,
    3.6373853495, -1.6467100004, -1.0946783260,
    3.6474846664, -0.1909460021, 1.9683378620
  };

  Data data = init_data(natoms, attyp, coord, charge, uhf);

  double* gradient = NULL;
  double energy = calc_energy_and_gradient(&data, &gradient);

  printf("Energy: %8.6f Eh\n", energy);
  if (gradient != NULL) {
    printf("Gradient vector per atom (Hartree/Bohr):\n");
    for (int i = 0; i < data.natoms; ++i) {
      printf("Atom %d (%2d): %12.6f %12.6f %12.6f\n",
            i+1,
            data.attyp[i],
            gradient[3*i + 0],
            gradient[3*i + 1],
            gradient[3*i + 2]);
    }


    free(gradient);
  }

  return 0;
}
