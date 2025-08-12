#include "convert.h"
#define TOL 1e-5

bool compare_xyz(const char* filename, int natoms,
                 char** symbols, double* xyz_data) {
  FILE* f = fopen(filename, "r");
  if (!f) {
    perror("Failed to open expected.xyz");
    return false;
  }

  char line[256];
  fgets(line, sizeof(line), f); // Skip first line
  fgets(line, sizeof(line), f); // Skip second line

  for (int i = 0; i < natoms; ++i) {
    char sym[8];
    double x, y, z;
    if (!fgets(line, sizeof(line), f)) {
      fprintf(stderr, "Expected file has fewer atoms.\n");
      fclose(f);
      return false;
    }
    sscanf(line, "%s %lf %lf %lf", sym, &x, &y, &z);

    if (strcmp(sym, symbols[i]) != 0) {
      fprintf(stderr, "Symbol mismatch at atom %d: %s vs %s\n", i, symbols[i], sym);
      fclose(f);
      return false;
    }

    for (int j = 0; j < 3; ++j) {
      double diff = fabs(xyz_data[3 * i + j] - ((j == 0) ? x : (j == 1) ? y : z));
      if (diff > TOL) {
        fprintf(stderr, "Coord mismatch at atom %d, coord %d: got %.6f, expected %.6f\n",
                i, j, xyz_data[3 * i + j], (j == 0) ? x : (j == 1) ? y : z);
        fclose(f);
        return false;
      }
    }
  }

  fclose(f);
  return true;
}

int main() {
  ZMatrix* zm = read_zmatrix("input.zmat");
  if (!zm) {
    fprintf(stderr, "Could not open file.\n");
    return 1;
  }

  int natoms = zm->natoms;

  // Allocate and fill geo[3 x natoms]
  double* geo = calloc(3 * natoms, sizeof(double));
  for (int i = 0; i < natoms; i++) {
    geo[0 * natoms + i] = zm->bond_len[i];
    geo[1 * natoms + i] = zm->bond_ang[i];
    geo[2 * natoms + i] = zm->dihed_ang[i];
  }

  // Allocate coord
  double* coord = calloc(3 * natoms, sizeof(double));
  if (!coord) {
    perror("calloc");
    free(geo);
    free_zmatrix(zm);
    return 1;
  }

  // Run geometry construction
  bool fail = gmetry(natoms, geo, zm->na, zm->nb, zm->nc, coord);
  if (fail) {
    fprintf(stderr, "Geometry construction failed.\n");
    free(coord);
    free(geo);
    free_zmatrix(zm);
    return 1;
  }

  // Created by: python3 zmattoxyz.py input.zmat
  if (!compare_xyz("expected.xyz", zm->natoms, zm->symbols, coord)) {
    fprintf(stderr, "Mismatch found in XYZ data\n");
  } else {
    printf("XYZ matches expected output!\n");
  }

  free(coord);
  free(geo);
  free_zmatrix(zm);
  return 0;
}
