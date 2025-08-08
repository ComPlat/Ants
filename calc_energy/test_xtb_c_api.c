#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "xtb.h"
static double POS_INF = 1.0 /0.0;

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

double calc_energy(Data* data) {
  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol = xtb_newMolecule(
    env, &(data->natoms), data->attyp, data->coord, &(data->charge), &(data->uhf), NULL, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return POS_INF;
  }

  xtb_setVerbosity(env,  XTB_VERBOSITY_MUTED);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return POS_INF;
  }

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return POS_INF;
  }

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return POS_INF;
  }

  double energy = 0.0;
  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return POS_INF;
  }

  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);

  return energy;
}

double* calc_gradient(Data* data) {
  double* gradient = NULL;
  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol = xtb_newMolecule(
    env, &(data->natoms), data->attyp, data->coord, &(data->charge), &(data->uhf), NULL, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return gradient;
  }

  xtb_setVerbosity(env,  XTB_VERBOSITY_MUTED);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return gradient;
  }

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return gradient;
  }

  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return gradient;
  }

  gradient = (double*) malloc(3 * data->natoms * sizeof(double));
  xtb_getGradient(env, res, gradient);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return gradient;
  }

  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);

  return gradient;
}

void backtrack(double* x, double* grad, double lr, double beta) {
  while(true) {

  }
}

double gradient_decent(double lr, double tol, int max_iter, double h) {

}

int main (int argc, char** argv) {

  int    const natoms = 18;
  int    const attyp[18] = {6, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1, 1};
  double const charge = 0.0;
  int    const uhf = 0;
  double const coord[3*18] = {
    0.0000000000, 0.0000000000, 0.0000000000,
    1.5517640000, 0.0000000000, 0.0000000000,
    -0.6214131653, 1.4212710069, 0.0000000000,
    -0.3202344715, -0.5777128156, 0.8955541160,
    -0.4020743303, -0.5228367172, -0.8938873510,
    0.1841653268, 2.3326703516, 0.9214999124,
    -1.6903042651, 1.3588550321, 0.2933374010,
    -0.5790544891, 1.8560721600, -1.0231949363,
    1.6441078582, 2.4786929260, 0.4164114954,
    0.1764956279, 1.8444387137, 1.9208853288,
    -0.2736817304, 3.3394546876, 1.0355205414,
    2.0956513638, 1.3275671641, -0.5201999921,
    2.2911111163, 2.5382257523, 1.3198049476,
    1.7800601271, 3.4231195994, -0.1525003770,
    1.6822536598, 1.4825678003, -1.5413943753,
    3.2027057981, 1.3293805938, -0.6032940337,
    1.9248779896, -0.8723317874, -0.5797090083,
    1.9302038454, -0.1010743891, 1.0410350959
  };
  Data data = init_data(natoms, attyp, coord, charge, uhf);

  double energy = calc_energy(&data);
  printf("%8.3f \n", energy);
  double* gradient = calc_gradient(&data);
  printf("%8.3f", gradient[0]);

  free(gradient);
  return 0;
}
