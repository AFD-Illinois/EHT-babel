#ifndef COORDINATES_H
#define COORDINATES_H

#include <math.h>
#include <string.h>

typedef int coords;
#define KORAL_BLCOORDS (1)
#define KSCOORDS       (2)
#define MKS3COORDS     (3)

coords get_coords_from_string(char *s);

// aliases
double Tan(double a);
double Power(double a, double b);
double Cos(double a);
double Csc(double a);
double Cot(double a);
double Sec(double a);

// metric stuff
extern double a;
extern double R0, H0, MY1, MY2, MP0;

// borrowed from koral_lite:metric.c
void coco_MKS32KS(double xMKS3[4], double xKS[4]);
void coco_KS2MKS3(double xKS[4], double xMKS3[4]);
void coco_BL2KS(double xBL[4], double xKS[4]);
void coco_KS2BL(double xKS[4], double xBL[4]);

void coco_EKS2KS(double xEKS[4], double xKS[4]);
void coco_KS2EKS(double xKS[4], double xEKS[4]);

void dxdx_KS2EKS(double xKS[4], double dxdxKS2EKS[4][4]);
void dxdx_EKS2KS(double xEKS[4], double dxdxEKS2KS[4][4]);

void dxdx_BL2KS(double xBL[4], double dxdxBL2KS[4][4]);
void dxdx_KS2BL(double xKS[4], double dxdxKS2BL[4][4]);

void dxdx_KS2MKS3(double xKS[4], double dxdxKS2MKS3[4][4]);
void dxdx_MKS32KS(double xMKS3[4], double dxdxMKS32KS[4][4]);

#endif // COORDINATES_H