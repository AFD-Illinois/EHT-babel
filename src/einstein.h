#ifndef EINSTEIN_H
#define EINSTEIN_H

// implements contractions
double a_a_(double u[4], double v[4]);
void ab_b_a(double u[4][4], double v[4], double w[4]);
void ab_a_b(double u[4][4], double v[4], double w[4]);

// potentially misleading name
void ik_jl_kl_ij(double T1[4][4], double T2[4][4], double A[4][4]);

#endif // EINSTEIN_H