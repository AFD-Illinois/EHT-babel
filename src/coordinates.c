
#include "coordinates.h"

double a=0.;
double R0=0., H0=0., MY1=0., MY2=0., MP0=0.;

double Pi = M_PI;

double Tan(double a) { return tan(a); }
double Power(double a, double b) { return pow(a,b); }
double Cos(double a) { return cos(a); }
double Csc(double a) { return 1/sin(a); }
double Cot(double a) { return 1/tan(a); }
double Sec(double a) { return 1/cos(a); }

coords get_coords_from_string(char *s) {
  // parse input string s to associated value

  if (strcmp(s, "MKS3COORDS")==0) return MKS3COORDS;
  else if (strcmp(s, "KSCOORDS")==0) return KSCOORDS;
  else if (strcmp(s, "BLCOORDS")==0) return KORAL_BLCOORDS;
  else if (strcmp(s, "BLCOORDS//KERRCOORDS")==0) return KORAL_BLCOORDS;
  else if (strcmp(s, "KERRCOORDS")==0) return KORAL_BLCOORDS;

  return 0;
}

void coco_MKS32KS(double xMKS3[4], double xKS[4]) {
  // translates from xMKS3 -> xKS

  double x0 = xMKS3[0];
  double x1 = xMKS3[1];
  double x2 = xMKS3[2];
  double x3 = xMKS3[3];

  double KSx0,KSx1,KSx2,KSx3;
  KSx0 = x0;
  KSx1 = exp(x1) + R0;
  KSx2 = (Pi*(1. + 1./tan((H0*Pi)/2.)*tan(H0*Pi*(-0.5 + (MY1 + (pow(2,MP0)*(-MY1 + MY2))/pow(exp(x1) + R0,MP0))*(1 - 2*x2) + x2))))/2.;
  KSx3 = x3;

  xKS[0] = KSx0;
  xKS[1] = KSx1;
  xKS[2] = KSx2;
  xKS[3] = KSx3;
}

void coco_KS2MKS3(double xKS[4], double xMKS3[4]) {
  // transform from xKS -> xMKS3
  double KSx0 = xKS[0];
  double KSx1 = xKS[1];
  double KSx2 = xKS[2];
  double KSx3 = xKS[3];

  double x0,x1,x2,x3;

  x0 = KSx0;  
  x1 = log(KSx1 - R0);
  x2 = (-(H0*pow(KSx1,MP0)*Pi) - pow(2,1 + MP0)*H0*MY1*Pi + 2*H0*pow(KSx1,MP0)*MY1*Pi + pow(2,1 + MP0)*H0*MY2*Pi + 2*pow(KSx1,MP0)*atan(((-2*KSx2 + Pi)*tan((H0*Pi)/2.))/Pi))/(2.*H0*(-pow(KSx1,MP0) - pow(2,1 + MP0)*MY1 + 2*pow(KSx1,MP0)*MY1 + pow(2,1 + MP0)*MY2)*Pi);
  x3 = KSx3;

  xMKS3[0] = x0;
  xMKS3[1] = x1;
  xMKS3[2] = x2;
  xMKS3[3] = x3;
}

void coco_BL2KS(double xBL[4], double xKS[4]) {
  // transform from xBL -> xKS
  double r = xBL[1];
  double delta = r*r - 2.*r + a*a;
  double sqrta = sqrt(1.-a*a);

  xKS[0] = xBL[0]+2./sqrta*atanh(sqrta/(1.-r))+log(delta);
  xKS[1] = xBL[1];
  xKS[2] = xBL[2];
  xKS[3] = xBL[3]; // sic KORAL
}

void coco_KS2BL(double xKS[4], double xBL[4]) {
  // transform from xKS -> xBL

  double r = xKS[1];
  double delta = r*r - 2.*r + a*a;
  double sqrta = sqrt(1.-a*a);

  xBL[0] = xKS[0]-2./sqrta*atanh(sqrta/(1.-r))-log(delta);
  xBL[1] = xKS[1];
  xBL[2] = xKS[2];
  xBL[3] = xKS[3]; // sic KORAL
}

void coco_EKS2KS(double xEKS[4], double xKS[4]) {
  // transform from xKS -> xBL

  xKS[0] = xEKS[0];
  xKS[1] = exp(xEKS[1]);
  xKS[2] = M_PI * xEKS[2];
  xKS[3] = xEKS[3];
}

void coco_KS2EKS(double xKS[4], double xEKS[4]) {
  // transform from xEKS -> xKS

  xEKS[0] = xKS[0];
  xEKS[1] = log(xKS[1]);
  xEKS[2] = xKS[2] / M_PI;
  xEKS[3] = xKS[3];
}


void dxdx_BL2KS(double xBL[4], double dxdxBL2KS[4][4]) {
  // get the Jacobian from BL -> KS at xBL (in BL coordinates)
  double t = xBL[0];
  double r = xBL[1];
  double th = xBL[2];
  double ph = xBL[3];

  double delta = r*r - 2.*r + a*a;

  for(int i=0; i<4; ++i) {
    for(int j=0; j<4; ++j) {
      dxdxBL2KS[i][j] = i==j ? 1 : 0;
    }
  }
  
  dxdxBL2KS[0][1] = 2.*r/delta;
  dxdxBL2KS[3][1] = a/delta;
}

void dxdx_KS2BL(double xKS[4], double dxdxKS2BL[4][4]) {
  // get the Jacobian from KS -> BL at xKS (in KS coordinates)
  double t = xKS[0];
  double r = xKS[1];
  double th = xKS[2];
  double ph = xKS[3];

  double delta = r*r - 2.*r + a*a;

  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      dxdxKS2BL[i][j] = i==j ? 1 : 0;
    }
  }
  
  dxdxKS2BL[0][1] = -2.*r/delta;
  dxdxKS2BL[3][1] = -a/delta;
}

void dxdx_KS2EKS(double xKS[4], double dxdxKS2EKS[4][4]) {
  // get the Jacobian from KS -> EKS at xKS (in KS coordinates)
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      dxdxKS2EKS[i][j] = i==j ? 1 : 0; 
    }
  }

  dxdxKS2EKS[1][1] = 1. / xKS[1];
  dxdxKS2EKS[2][2] = 1. / M_PI;
}

void dxdx_EKS2KS(double xEKS[4], double dxdxEKS2KS[4][4]) {
  // get the Jacobian from EKS -> KS at xEKS (in EKS coordinates)
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      dxdxEKS2KS[i][j] = i==j ? 1 : 0; 
    }
  }

  dxdxEKS2KS[1][1] = exp(xEKS[1]);
  dxdxEKS2KS[2][2] = M_PI;
}

void dxdx_KS2MKS3(double xKS[4], double dxdxKS2MKS3[4][4]) {
  // get the Jacobian from KS -> MKS3 at xKS (in KS coordinates)
  double KSx0 = xKS[0];
  double KSx1 = xKS[1];
  double KSx2 = xKS[2];
  double KSx3 = xKS[3];

  for (int i=0; i<4; ++i) {
    for(int j=0; j<4; ++j) {
      dxdxKS2MKS3[i][j] = i==j ? 1 : 0;
    }
  }
  
  dxdxKS2MKS3[1][1]= 1./(KSx1 - R0);
  dxdxKS2MKS3[2][1]= -((Power(2,1 + MP0)*Power(KSx1,-1 + MP0)*MP0*(MY1 - MY2)*atan(((-2*KSx2 + Pi)*Tan((H0*Pi)/2.))/Pi))/(H0*Power(Power(KSx1,MP0)*(1 - 2*MY1) + Power(2,1 + MP0)*(MY1 - MY2),2)*Pi));
  dxdxKS2MKS3[2][2]= (-2*Power(KSx1,MP0)*Tan((H0*Pi)/2.))/(H0*(Power(KSx1,MP0)*(-1 + 2*MY1) + Power(2,1 + MP0)*(-MY1 + MY2))*Power(Pi,2)*(1 + (Power(-2*KSx2 + Pi,2)*Power(Tan((H0*Pi)/2.),2))/Power(Pi,2)));
}


void dxdx_MKS32KS(double xMKS3[4], double dxdxMKS32KS[4][4]) {
  // get the Jacobian from MKS3 -> KS at xMKS3 (in MKS3 coordinates)
  double x0 = xMKS3[0];
  double x1 = xMKS3[1];
  double x2 = xMKS3[2];
  double x3 = xMKS3[3];

  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      dxdxMKS32KS[i][j] = i==j ? 1 : 0;
    }
  }
  
  dxdxMKS32KS[1][1]= exp(x1);
  dxdxMKS32KS[2][1]= -(Power(2,-1 + MP0)*exp(x1)*H0*MP0*(MY1 - MY2)*Power(Pi,2)*Power(exp(x1) + R0,-1 - MP0)*(-1 + 2*x2)*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2));
  dxdxMKS32KS[2][2]= (H0*Power(Pi,2)*(1 - 2*(MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0)))*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2))/2.;
}

