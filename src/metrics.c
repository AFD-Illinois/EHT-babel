#include "metrics.h"

void calc_gcon_mks3(double *xMKS3, double gcon[4][4]) {
  // returns contravariant metric at point xMKS3 (in MKS3 coordinates) with components in MKS3 basis

  // first get gcon in KS basis at xMKS3
  double xKS[4], gconKS[4][4] = { 0 };
  coco_MKS32KS(xMKS3, xKS);
  calc_gcon_ks(xKS, gconKS);

  // then transform to MKS3 basis
  double dxdx[4][4] = { 0 };
  dxdx_KS2MKS3(xKS, dxdx);
  ik_jl_kl_ij(gconKS, gcon, dxdx);
}

void calc_gcon_eks(double *xEKS, double gcon[4][4]) {
  // returns contravariant metric at point xEKS (in EKS coordinates) with components in EKS basis

  // first get gcon in KS basis at xMKS3
  double xKS[4], gconKS[4][4] = { 0 };
  coco_EKS2KS(xEKS, xKS);
  calc_gcon_ks(xKS, gconKS);

  // then transform to MKS3 basis
  double dxdx[4][4] = { 0 };
  dxdx_KS2EKS(xKS, dxdx);
  ik_jl_kl_ij(gconKS, gcon, dxdx);
}


void calc_gcon_bl(double *xBL, double gcon[4][4]) {
  // return contravariant metric at point xBL (in BL coordinates) with components in BL basis

  double x0 = xBL[0];
  double x1 = xBL[1];
  double x2 = xBL[2];
  double x3 = xBL[3];

  gcon[0][0]= -((Power(a,4) + 2*Power(x1,4) + Power(a,2)*x1*(2 + 3*x1) + Power(a,2)*(Power(a,2) + (-2 + x1)*x1)*Cos(2*x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))));
  gcon[0][1]= 0;
  gcon[0][2]= 0;
  gcon[0][3]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)));
  gcon[1][0]= 0;
  gcon[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[1][2]= 0;
  gcon[1][3]= 0;
  gcon[2][0]= 0;
  gcon[2][1]= 0;
  gcon[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[2][3]= 0;
  gcon[3][0]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)));
  gcon[3][1]= 0;
  gcon[3][2]= 0;
  gcon[3][3]= (2*((-2 + x1)*x1 + Power(a,2)*Power(Cos(x2),2))*Power(Csc(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)));
}

void calc_gcon_ks(double *xKS, double gcon[4][4]) {
  // return contravariant metric at point xKS (in KS coordinates) with components in KS basis

  double x0 = xKS[0];
  double x1 = xKS[1];
  double x2 = xKS[2];
  double x3 = xKS[3];

  gcon[0][0] = -((x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)));
  gcon[0][1] = (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[0][2] = 0;
  gcon[0][3] = 0;
  gcon[1][0] = (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[1][1] = (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[1][2] = 0;
  gcon[1][3] = a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[2][0] = 0;
  gcon[2][1] = 0;
  gcon[2][2] = 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[2][3] = 0;
  gcon[3][0] = 0;
  gcon[3][1] = a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
  gcon[3][2] = 0;
  gcon[3][3] = Power(Csc(x2),2)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2));
}


// taken from koral_lite:misc.c
int invert_44(double a[][4], double ia[][4]) {

  double mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  double tmp[12];
  double src[16];
  double det;
  
  // transpose matrix
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  
  // calculate pairs for first 8 elements (cofactors)
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  
  // calculate first 8 elements (cofactors)
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];

  // calculate pairs for second 8 elements (cofactors)
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];

  // calculate second 8 elements (cofactors)
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];

  // calculate determinant
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

  // calculate matrix inverse
  det = 1/det; 

  if(isnan(det))
    //my_err("det in inverse 4x4 zero\n");
    return -1;

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++){
      ia[i][j]= dst[i*4+j];
      if(isnan(ia[i][j])) return -1;
    }
  }

  return 0;
}
