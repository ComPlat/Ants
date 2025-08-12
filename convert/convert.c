#include "convert.h"

typedef struct {
  double* data;
  int nrow;
  int ncol;
} Array;

double get(Array* array, int row, int col) {
  return array->data[(row * array->ncol) + col];
}
void set(Array* array, int row, int col, double val) {
  array->data[(row * array->ncol) + col] = val;
}

bool gmetry(int natoms, double* geo_data,
            int* na, int* nb, int* nc,
            double* xyz) {
  bool fail = false;

  double* coord_data = (double*) calloc(3 * natoms, sizeof(double));
  if (!coord_data) return true;

  Array coord = { .data = coord_data, .nrow = 3, .ncol = natoms };
  Array geo   = { .data = geo_data,   .nrow = 3, .ncol = natoms };

  // Initialize first two atoms
  set(&coord, 0, 0, 0.0);
  set(&coord, 1, 0, 0.0);
  set(&coord, 2, 0, 0.0);

  set(&coord, 0, 1, get(&geo, 0, 1));
  set(&coord, 1, 1, 0.0);
  set(&coord, 2, 1, 0.0);

  if (natoms == 2) {
    memcpy(xyz, coord_data, sizeof(double) * 3 * natoms);
    free(coord_data);
    return false;
  }

  // Atom 3
  double angle = get(&geo, 1, 2);
  double ccos = cos(angle);

  if (na[2] == 0)
    set(&coord, 0, 2, get(&coord, 0, 0) + get(&geo, 0, 2) * ccos);
  else
    set(&coord, 0, 2, get(&coord, 0, 1) - get(&geo, 0, 2) * ccos);

  set(&coord, 1, 2, get(&geo, 0, 2) * sin(angle));
  set(&coord, 2, 2, 0.0);

  // Atoms >= 4
  for (int i = 3; i < natoms; i++) {
    double cosa = cos(get(&geo, 1, i));
    int mb = nb[i];
    int mc = na[i];

    double xb = get(&coord, 0, mb) - get(&coord, 0, mc);
    double yb = get(&coord, 1, mb) - get(&coord, 1, mc);
    double zb = get(&coord, 2, mb) - get(&coord, 2, mc);

    double rbc = 1.0 / sqrt(xb * xb + yb * yb + zb * zb);

    if (fabs(cosa) >= 0.9999999999) {
      double rbc2 = get(&geo, 0, i) * rbc * cosa;
      set(&coord, 0, i, get(&coord, 0, mc) + xb * rbc2);
      set(&coord, 1, i, get(&coord, 1, mc) + yb * rbc2);
      set(&coord, 2, i, get(&coord, 2, mc) + zb * rbc2);
      continue;
    }

    int ma = nc[i];

    double xa = get(&coord, 0, ma) - get(&coord, 0, mc);
    double ya = get(&coord, 1, ma) - get(&coord, 1, mc);
    double za = get(&coord, 2, ma) - get(&coord, 2, mc);

    double xyb = sqrt(xb * xb + yb * yb);
    int k = -1;
    if (xyb <= 0.1) {
      double xpa = za; za = -xa; xa = xpa;
      double xpb = zb; zb = -xb; xb = xpb;
      xyb = sqrt(xb * xb + yb * yb);
      k = 1;
    }

    double costh = xb / xyb;
    double sinth = yb / xyb;

    double xpa = xa * costh + ya * sinth;
    double ypa = ya * costh - xa * sinth;

    double sinph = zb * rbc;
    double cosph = sqrt(fabs(1.0 - sinph * sinph));

    double xqa = xpa * cosph + za * sinph;
    double zqa = za * cosph - xpa * sinph;

    double yza = sqrt(ypa * ypa + zqa * zqa);
    if (yza < 1e-4) {
      fail = true;
      break;
    }

    double coskh = ypa / yza;
    double sinkh = zqa / yza;

    double sina = sin(get(&geo, 1, i));
    double sind = -sin(get(&geo, 2, i));
    double cosd = cos(get(&geo, 2, i));

    double xd = get(&geo, 0, i) * cosa;
    double yd = get(&geo, 0, i) * sina * cosd;
    double zd = get(&geo, 0, i) * sina * sind;

    double ypd = yd * coskh - zd * sinkh;
    double zpd = zd * coskh + yd * sinkh;
    double xpd = xd * cosph - zpd * sinph;
    double zqd = zpd * cosph + xd * sinph;
    double xqd = xpd * costh - ypd * sinth;
    double yqd = ypd * costh + xpd * sinth;

    if (k < 1) {
      set(&coord, 0, i, xqd + get(&coord, 0, mc));
      set(&coord, 1, i, yqd + get(&coord, 1, mc));
      set(&coord, 2, i, zqd + get(&coord, 2, mc));
    } else {
      set(&coord, 0, i, -zqd + get(&coord, 0, mc));
      set(&coord, 1, i, yqd + get(&coord, 1, mc));
      set(&coord, 2, i, xqd + get(&coord, 2, mc));
    }
  }

  for (int i = 0; i < natoms; i++) {
    xyz[i * 3 + 0] = get(&coord, 0, i);
    xyz[i * 3 + 1] = get(&coord, 1, i);
    xyz[i * 3 + 2] = get(&coord, 2, i);
  }
  free(coord_data);
  return fail;
}

int symbol_to_atomic_number(const char* symbol) {
  if (strcmp(symbol, "H") == 0) return 1;
  if (strcmp(symbol, "C") == 0) return 6;
  if (strcmp(symbol, "N") == 0) return 7;
  if (strcmp(symbol, "O") == 0) return 8;
  if (strcmp(symbol, "F") == 0) return 9;
  if (strcmp(symbol, "Cl") == 0) return 17;
  if (strcmp(symbol, "Br") == 0) return 35;
  if (strcmp(symbol, "I") == 0) return 53;
  // Add more as needed
  return -1; // Unknown
}

void free_zmatrix(ZMatrix* zm) {
  if (!zm) return;
  for (int i = 0; i < zm->natoms; i++) {
    free(zm->symbols[i]);
  }
  free(zm->symbols);
  free(zm->na);
  free(zm->nb);
  free(zm->nc);
  free(zm->bond_len);
  free(zm->bond_ang);
  free(zm->dihed_ang);
  free(zm->atomic_numbers);
  free(zm);
}

ZMatrix* read_zmatrix(const char* filename) {
  FILE* f = fopen(filename, "r");
  if (!f) return NULL;

  ZMatrix* zm = calloc(1, sizeof(ZMatrix));
  int capacity = 8;

  zm->symbols    = malloc(capacity * sizeof(char*));
  zm->na         = malloc(capacity * sizeof(int));
  zm->nb         = malloc(capacity * sizeof(int));
  zm->nc         = malloc(capacity * sizeof(int));
  zm->bond_len   = malloc(capacity * sizeof(double));
  zm->bond_ang   = malloc(capacity * sizeof(double));
  zm->dihed_ang  = malloc(capacity * sizeof(double));
  zm->atomic_numbers = malloc(capacity * sizeof(int));

  char* line = NULL;
  size_t len = 0;

  while (getline(&line, &len, f) != -1) {
    // Skip empty lines
    char* trimmed = strtok(line, "\n");
    if (!trimmed || strlen(trimmed) == 0) continue;

    // Resize arrays if needed
    if (zm->natoms >= capacity) {
      capacity *= 2;
      zm->symbols   = realloc(zm->symbols, capacity * sizeof(char*));
      zm->na        = realloc(zm->na, capacity * sizeof(int));
      zm->nb        = realloc(zm->nb, capacity * sizeof(int));
      zm->nc        = realloc(zm->nc, capacity * sizeof(int));
      zm->bond_len  = realloc(zm->bond_len, capacity * sizeof(double));
      zm->bond_ang  = realloc(zm->bond_ang, capacity * sizeof(double));
      zm->dihed_ang = realloc(zm->dihed_ang, capacity * sizeof(double));
      zm->atomic_numbers = realloc(zm->atomic_numbers, capacity * sizeof(int));
    }

    char atom[16];
    int na = 0, nb = 0, nc = 0;
    double bl = 0.0, ba = 0.0, da = 0.0;

    int i = zm->natoms;
    if (i == 0) {
      sscanf(trimmed, "%s", atom);
    } else if (i == 1) {
      sscanf(trimmed, "%s %d %lf", atom, &na, &bl);
    } else if (i == 2) {
      sscanf(trimmed, "%s %d %lf %d %lf", atom, &na, &bl, &nb, &ba);
      ba *= M_PI / 180.0;
    } else {
      sscanf(trimmed, "%s %d %lf %d %lf %d %lf", atom, &na, &bl, &nb, &ba, &nc, &da);
      ba *= M_PI / 180.0;
      da *= M_PI / 180.0;
    }

    zm->symbols[i]    = strdup(atom);
    zm->atomic_numbers[i] = symbol_to_atomic_number(atom);
    zm->na[i]         = na > 0 ? na - 1 : 0;
    zm->nb[i]         = nb > 0 ? nb - 1 : 0;
    zm->nc[i]         = nc > 0 ? nc - 1 : 0;
    zm->bond_len[i]   = bl;
    zm->bond_ang[i]   = ba;
    zm->dihed_ang[i]  = da;

    zm->natoms++;
  }

  free(line);
  fclose(f);
  return zm;
}
