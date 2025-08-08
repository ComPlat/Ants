#include "calc_energy.h"

static double POS_INF = 1.0 /0.0;

void init_settings_loss_fct(settings_loss_fct* slf, const int natoms, const int* attyp,
                            int* na, int* nb, int* nc,
                            double electronic_temperature,
                            double accuracy, int max_iter, bool verbose) {
  slf->natoms = natoms;
  slf->attyp = (int*) malloc(natoms * sizeof(int));
  memcpy(slf->attyp, attyp, natoms * sizeof(int));
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
  slf->coord = (double*) calloc(3 * slf->natoms, sizeof(double));

  slf->na = malloc(natoms * sizeof(int));
  slf->nb = malloc(natoms * sizeof(int));
  slf->nc = malloc(natoms * sizeof(int));

  memcpy(slf->na, na, natoms * sizeof(int));
  memcpy(slf->nb, nb, natoms * sizeof(int));
  memcpy(slf->nc, nc, natoms * sizeof(int));
}


void destroy_settings_loss_fct(settings_loss_fct* slf) {
  xtb_delete(slf->res);
  xtb_delete(slf->calc);
  xtb_delete(slf->mol);
  xtb_delete(slf->env);
  free(slf->attyp); slf->attyp = NULL;
  free(slf->coord); slf->coord = NULL;
  free(slf->na); slf->na = NULL;
  free(slf->nb); slf->nb = NULL;
  free(slf->nc); slf->nc = NULL;
}

double calc_energy(double* coord, void* data) {
  settings_loss_fct* slf = (settings_loss_fct*)data;
  bool fail = gmetry(slf->natoms, coord, slf->na, slf->nb, slf->nc, slf->coord);
  for (int i = 0; i < slf->natoms; i++) {
    printf("%d ", slf->attyp[i]); // correct
  }
  printf("\n");
  for (int i = 0; i < (3*slf->natoms); i++) {
    printf("%8.6f ", slf->coord[i]); // correct
    if ((i + 1) % 3 == 0) printf("\n");
  }
  printf("\n");
  if (fail) {
    return POS_INF;
  }

  double const charge = 0.0; // cyclo hexan
  int    const uhf = 0; // what is this???
  if (slf->mol == NULL) {
    slf->mol = xtb_newMolecule(
      slf->env,
      &slf->natoms,
      slf->attyp,
      slf->coord,
      &charge,
      // NULL, // const double* charge in e
      &uhf,
      // NULL, // const int* uhf
      NULL, // const double* lattice[3][3]
      NULL // const bool* periodic [3]
    );
  } else {
    xtb_updateMolecule(slf->env, slf->mol, slf->coord, NULL);
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


double run_xtb_energy(double* coord, void* data) {
  ZMatrix* zm = (ZMatrix*)data;
  double* geo = calloc(3 * zm->natoms, sizeof(double));
  bool fail = gmetry(zm->natoms, coord, zm->na, zm->nb, zm->nc, geo);
  if (fail) {
    return POS_INF;
  }

  double energy = 0.0;

  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();

  const int uhf = 0;
  const double charge = 0.0;
  xtb_TMolecule mol = xtb_newMolecule(env, &zm->natoms, zm->atomic_numbers, geo, &charge, &uhf, NULL, NULL);

  xtb_loadGFN2xTB(env, mol, calc, NULL);
  xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);

  xtb_singlepoint(env, mol, calc, res);
  xtb_getEnergy(env, res, &energy);

  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);
  free(geo);

  return energy;
}
