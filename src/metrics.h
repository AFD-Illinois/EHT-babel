#ifndef METRICS_H
#define METRICS_H

#include <math.h>

#include "einstein.h"
#include "coordinates.h" // for metric parameters

// numerical contravariant metric calculations
void calc_gcon_mks3(double *xMKS3, double gcon[4][4]);
void calc_gcon_eks(double *xEKS, double gcon[4][4]);

// analytic contravariant metric calculations
void calc_gcon_bl(double *bl, double gcon[4][4]);
void calc_gcon_ks(double *ks, double gcon[4][4]);

// linear algebra
int invert_44(double a[][4], double ia[][4]);

#endif // METRICS_H