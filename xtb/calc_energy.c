#include "calc_energy.h"

static double POS_INF = 1.0 /0.0;

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

double calc_energy(const double* coord, void* data) {

  settings_loss_fct* slf = (settings_loss_fct*)data;
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
