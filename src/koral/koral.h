#ifndef KORAL_H
#define KORAL_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../utils.h"
#include "../coordinates.h"
#include "../metrics.h"
#include "../einstein.h"

#define P4VEC(S,X) fprintf(stderr, "%s: %g %g %g %g\n", S, X[0], X[1], X[2], X[3]);

// these could change
#define KORAL_XI_OFFSET (2)
#define MU_GAS (0.5)  // see below

// the full details of MU_GAS are given here and come from choices.h
// #define MU_GAS (1./(0.5 + 1.5*HFRAC + 0.25*HEFRAC + A_INV_MEAN*MFRAC))
// #define HFRAC 1.
// #define MFRAC (1.-HFRAC-HEFRAC)
// #define HEFRAC 0.
// #define A_INV_MEAN 0.0570490404519

// these shouldn't change  
#define K_BOLTZ_CGS (1.3806488e-16)
#define M_PROTON_CGS (1.67262158e-24)
#define M_ELECTRON_CGS (9.1093826e-28)
#define CL_CGS (2.99792458e10)

// KORAL-internal units. used for conversion between CODES & CGS
// the following list has been abbreviated to remove tilde units
#define KORAL_CCC (2.9979246e10)
#define KORAL_GGG (6.674e-8)
#define KORAL_MSUNCM (147700.)

// KORAL outputs different file formats
typedef int koral_dtype;
#define KORAL_GRMHD (1)
#define KORAL_RADGRMHD (2)

typedef struct {
  coords zone_coordinates;
  coords file_coordinates;
  double a;
  double R0;
  double H0;
  double MY1;
  double MY2;
  double MP0;
  int nx;
  int ny;
  int nz;
  int np;

  double t;
  double mass;
  double rhoCGS2GU;
  double uintCGS2GU;
  double endenCGS2GU;

  double startx[4];
  double dx[4];

  koral_dtype dtype;
  double DTOUT1;
  double gam;
  double gam_e;
  double gam_p;
} geom_koral;

int koral_init(char *fname, geom_koral *geom);
int koral_read(char *fname, double ****data, geom_koral *geom);
int koral_translate(double ****rawdata, double *prims, geom_koral *geom);

int koral_init_RadGRMHD(char *fname, geom_koral *geom);
int koral_init_GRMHD(char *fname, geom_koral *geom);

void set_u0(double ucon[4], double gcon[4][4], double gcov[4][4]);

#endif // KORAL_H
