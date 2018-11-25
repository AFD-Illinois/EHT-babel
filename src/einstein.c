#include "einstein.h"

double a_a_(double u[4], double v[4]) {
  // u^a v_a =: retval
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3];
}

void ab_b_a(double u[4][4], double v[4], double w[4]) {
  // u^{ab} v_{b} =: w^a
  w[0] = 0.; 
  w[1] = 0.;
  w[2] = 0.;
  w[3] = 0.;

  for (int i=0; i<4; ++i) {
    w[0] += u[0][i] * v[i];
    w[1] += u[1][i] * v[i];
    w[2] += u[2][i] * v[i];
    w[3] += u[3][i] * v[i];
  }
}

void ab_a_b(double u[4][4], double v[4], double w[4]) {
  // u^{ab} v_{a} =: w^b
  w[0] = 0.; 
  w[1] = 0.;
  w[2] = 0.;
  w[3] = 0.;

  for (int i=0; i<4; ++i) {
    w[0] += u[i][0] * v[i];
    w[1] += u[i][1] * v[i];
    w[2] += u[i][2] * v[i];
    w[3] += u[i][3] * v[i];
  }
}

// borrowed in part from koral:frames.c
void ik_jl_kl_ij(double T1[4][4], double T2[4][4], double A[4][4]) {
  // T2^{ij} = A^i_k A^j_l T1^{kl}
  int i,j,k,l;
  
  double Tt[4][4];

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      Tt[i][j]=T1[i][j];
    }
  }

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      T2[i][j]=0.;
      for(k=0;k<4;k++) {
        for(l=0;l<4;l++) {
          T2[i][j]+=A[i][k]*A[j][l]*Tt[k][l];
        }
      }
    }
  }
}

