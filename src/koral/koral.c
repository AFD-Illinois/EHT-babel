#include "koral.h"

int koral_init(char *fname, geom_koral *geom) {
  // swap off different define.h readers based on the type of run

  int rval = 0;

  if (geom->dtype == KORAL_GRMHD) {
    if ((rval = koral_init_GRMHD(fname, geom)) != 0) return rval;
  } else if (geom->dtype == KORAL_RADGRMHD) {
    if ((rval = koral_init_RadGRMHD(fname, geom)) != 0) return rval;
  } else {
    return -1;
  }

  // for now, these are constants
  geom->gam_e = 4./3;
  geom->gam_p = 5./3;

  // set conversion units
  double MASSCM = geom->mass * KORAL_MSUNCM;
  geom->rhoCGS2GU = (KORAL_GGG*MASSCM*MASSCM) / (KORAL_CCC*KORAL_CCC);
  double K_BOLTZ = K_BOLTZ_CGS * KORAL_GGG / KORAL_CCC / KORAL_CCC / KORAL_CCC / KORAL_CCC / MASSCM;
  double M_PROTON = M_PROTON_CGS * KORAL_GGG / (MASSCM * KORAL_CCC * KORAL_CCC);
  geom->uintCGS2GU = geom->rhoCGS2GU * K_BOLTZ / M_PROTON / MU_GAS;
  geom->endenCGS2GU = (KORAL_GGG * MASSCM * MASSCM) / (KORAL_CCC * KORAL_CCC * KORAL_CCC * KORAL_CCC);

  return 0;
}

int koral_init_RadGRMHD(char *fname, geom_koral *geom) {
  // rather than waste time parsing, assume the following will be
  // defined somewhere in the file
  // BHSPIN, MKSR0, MKSH0, MKSMY1, MKSMY2, MKSMP0, MYCOORDS,
  // DTOUT1, OUTCOORDS, GAMMA

  geom->np = 16; // 3 coordinates + 2+4+4 "prims" + Te,Tp + Gamma

  FILE *fp = NULL;
  char *line = NULL;
  char op[256], v1[256], v2[256];
  size_t len = 0;
  ssize_t read = 0;

  fp = fopen(fname, "r");
  if (fp == NULL) return -1;

  while ((read = getline(&line, &len, fp)) != -1) {
    char *trimmed = trimstr(line);
    if (strlen(trimmed) < 1) continue;
    if (trimmed[0] == '/') continue;
    sscanf(trimmed, "%s %s %s", op, v1, v2);
    if (strcmp(op,"#define")==0) {
      if (strcmp(v1,"BHSPIN")==0) {
        geom->a = atof(v2);
      } else if (strcmp(v1,"MASS")==0) {
        geom->mass = atof(v2);
      } else if (strcmp(v1,"MKSR0")==0) {
        geom->R0 = atof(v2);
      } else if (strcmp(v1,"MKSH0")==0) {
        geom->H0 = atof(v2);
      } else if (strcmp(v1,"MKSMY1")==0) {
        geom->MY1 = atof(v2);
      } else if (strcmp(v1,"MKSMY2")==0) {
        geom->MY2 = atof(v2);
      } else if (strcmp(v1,"MKSMP0")==0) {
        geom->MP0 = atof(v2);
      } else if (strcmp(v1,"TNX")==0) {
        geom->nx = atoi(v2) + KORAL_XI_OFFSET;
      } else if (strcmp(v1,"TNY")==0) {
        geom->ny = atoi(v2);
      } else if (strcmp(v1,"TNZ")==0) {
        geom->nz = atoi(v2);
      } else if (strcmp(v1,"MYCOORDS")==0) {
        geom->zone_coordinates = get_coords_from_string(v2);
      } else if (strcmp(v1,"DTOUT1")==0) {
        geom->DTOUT1 = atof(v2);
      } else if (strcmp(v1,"OUTCOORDS")==0) {
        geom->file_coordinates = get_coords_from_string(v2);
      } else if (strcmp(v1,"GAMMA")==0) {
        if (strcmp(v2,"(5./3.)")==0) {
          geom->gam = 5./3.;
        } else if (strcmp(v2,"(13./9.)")==0) {
          geom->gam = 13./9.;
        } else {
          fprintf(stderr, "encountered unknown fluid gamma: %s\n", v2);
        }
      }
    }
  }

  fclose(fp);

  // set metric parameters. these live in coordinates.c (extern'd in .h)
  a = geom->a;
  R0 = geom->R0;
  H0 = geom->H0;
  MY1 = geom->MY1;
  MY2 = geom->MY2;
  MP0 = geom->MP0;

  // other stuff (TODO)
  geom->t = 0.;

  return 0;
}

int koral_init_GRMHD(char *fname, geom_koral *geom) {
  // rather than waste time parsing, assume the following will be
  // defined somewhere in the file
  // BHSPIN, MKSR0, MKSH0, MKSMY1, MKSMY2, MKSMP0, MYCOORDS,
  // DTOUT1, OUTCOORDS, GAMMA

  geom->np = 13; // 3 coordinates + 2+4+4 "prims"

  FILE *fp = NULL;
  char *line = NULL;
  char op[256], v1[256], v2[256];
  size_t len = 0;
  ssize_t read = 0;

  fp = fopen(fname, "r");
  if (fp == NULL) return -1;

  // rudimentary parser
  strlist_el *defined = NULL;
  int doread = 1;
  int within = 0;

  // this parser is very broken but it works for the format of file
  // given. TODO should look at this more in the future...
  while ((read = getline(&line, &len, fp)) != -1) {
    char *trimmed = trimstr(line);
    if (strlen(trimmed) < 1) continue;
    if (trimmed[0] == '/') continue;
    sscanf(trimmed, "%s %s %s", op, v1, v2);
    if (doread && strcmp(op,"#define")==0) {
      defined = push_strlist(v1, defined);
      if (strcmp(v1,"BHSPIN")==0) {
        geom->a = atof(v2);
      } else if (strcmp(v1,"MASS")==0) {
        geom->mass = atof(v2);
      } else if (strcmp(v1,"MKSR0")==0) {
        geom->R0 = atof(v2);
      } else if (strcmp(v1,"MKSH0")==0) {
        geom->H0 = atof(v2);
      } else if (strcmp(v1,"MKSMY1")==0) {
        geom->MY1 = atof(v2);
      } else if (strcmp(v1,"MKSMY2")==0) {
        geom->MY2 = atof(v2);
      } else if (strcmp(v1,"MKSMP0")==0) {
        geom->MP0 = atof(v2);
      } else if (strcmp(v1,"TNX")==0) {
        geom->nx = atoi(v2) + KORAL_XI_OFFSET;
      } else if (strcmp(v1,"TNY")==0) {
        geom->ny = atoi(v2);
      } else if (strcmp(v1,"TNZ")==0) {
        geom->nz = atoi(v2);
      } else if (strcmp(v1,"MYCOORDS")==0) {
        geom->zone_coordinates = get_coords_from_string(v2);
      } else if (strcmp(v1,"DTOUT1")==0) {
        geom->DTOUT1 = atof(v2);
      } else if (strcmp(v1,"OUTCOORDS")==0) {
        geom->file_coordinates = get_coords_from_string(v2);
      } else if (strcmp(v1,"GAMMA")==0) {
        if (strcmp(v2,"(5./3.)")==0) {
          geom->gam = 5./3.;
        } else if (strcmp(v2,"(4./3.)")==0) {
          geom->gam = 4./3.;
        } else {
          fprintf(stderr, "encountered unknown fluid gamma: %s\n", v2);
        }
      }
    } else if (strcmp(op,"#ifdef")==0) {
      if (strlist_count(v1, defined) == 0) {
        doread = 0;
        within += 1;
      }
    } else if (strcmp(op,"#endif")==0) {
      if (within > 0) within -= 1;
      if (within == 0) doread = 1;
    }
  }

  fclose(fp);

  delete_strlist(defined);
  defined = NULL;
  
  // set metric parameters. these live in coordinates.c (extern'd in .h)
  a = geom->a;
  R0 = geom->R0;
  H0 = geom->H0;
  MY1 = geom->MY1;
  MY2 = geom->MY2;
  MP0 = geom->MP0;

  // TODO maybe parse filename?
  geom->t = 0.;

  return 0;
}

int koral_read(char *fname, double ****data, geom_koral *geom) {

  FILE *fp = NULL;
  char *line = NULL;
  char op[256], v1[256], v2[256];
  size_t len = 0;
  ssize_t read = 0;

  fp = fopen(fname, "r");
  if (fp == NULL) return -1;

  if (geom->dtype == KORAL_GRMHD) {

    int idum, xi, xj, xk;
    double r, h, p;
    double rho, tgas, u0, u1, u2, u3, bsq, b1, b2, b3;
    double fdum;

    while ( fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &xi, &xj, &xk, &r, &h, &p,
      &rho, &tgas, &u0, &u1, &u2, &u3, &fdum, &bsq, &b1, &b2, &b3, 
      &fdum, &fdum, &fdum, &fdum) > 0 ) {

      xi += KORAL_XI_OFFSET;

      data[xi][xj][xk][0] = r;
      data[xi][xj][xk][1] = h;
      data[xi][xj][xk][2] = p;

      data[xi][xj][xk][3] = rho;
      data[xi][xj][xk][4] = tgas;

      data[xi][xj][xk][5] = u0;
      data[xi][xj][xk][6] = u1;
      data[xi][xj][xk][7] = u2;
      data[xi][xj][xk][8] = u3;

      data[xi][xj][xk][9] = bsq;
      data[xi][xj][xk][10] = b1;
      data[xi][xj][xk][11] = b2;
      data[xi][xj][xk][12] = b3;

    }

  } else if (geom->dtype == KORAL_RADGRMHD) {
    // first line is header information. TODO maybe check this against loaded geom?
    fscanf(fp, "%*[^\n]\n", NULL);

    int idum, xi, xj, xk;
    double r, h, p;
    double rho, tgas, u0, u1, u2, u3, b0, b1, b2, b3, te, ti, gamma;
    double fdum;

    while ( fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
        &xi, &xj, &xk, &fdum, &fdum, &fdum, &r, &h, &p,
        &rho, &tgas, &u0, &u1, &u2, &u3, &b0, &b1, &b2, &b3,
        &fdum, &te, &ti, &gamma) > 0 ) {

      xi += KORAL_XI_OFFSET;

      data[xi][xj][xk][0] = r;
      data[xi][xj][xk][1] = h;
      data[xi][xj][xk][2] = p;

      data[xi][xj][xk][3] = rho;
      data[xi][xj][xk][4] = tgas;

      data[xi][xj][xk][5] = u0;
      data[xi][xj][xk][6] = u1;
      data[xi][xj][xk][7] = u2;
      data[xi][xj][xk][8] = u3;

      data[xi][xj][xk][9] = b0;
      data[xi][xj][xk][10] = b1;
      data[xi][xj][xk][11] = b2;
      data[xi][xj][xk][12] = b3;

      data[xi][xj][xk][13] = te;
      data[xi][xj][xk][14] = ti;

      data[xi][xj][xk][15] = gamma;

    }

  }

  // get logical geometry
  int n1 = geom->nx-1;
  int n2 = geom->ny-1;
  int n3 = geom->nz-1;

  // technically these are KORAL_BLCOORDS, but they're the same as KS where we care.
  double xleftKS[4] = { 0., data[0][0][0][0], data[0][0][0][1], data[0][0][0][2] };
  double xrightKS[4] = { 0., data[n1][n2][n3][0], data[n1][n2][n3][1], data[n1][n2][n3][2] };

  // get MKS3 equivalent
  double xleftMKS3[4] = { 0. };
  double xrightMKS3[4] = { 0. };
  coco_KS2MKS3(xleftKS, xleftMKS3);
  coco_KS2MKS3(xrightKS, xrightMKS3);

  geom->dx[1] = ( xrightMKS3[1] - xleftMKS3[1] ) / n1;
  geom->dx[2] = ( xrightMKS3[2] - xleftMKS3[2] ) / n2;
  geom->dx[3] = ( xrightMKS3[3] - xleftMKS3[3] ) / n3;

  geom->startx[1] = xleftMKS3[1] - geom->dx[1]/2.;
  geom->startx[2] = xleftMKS3[2] - geom->dx[2]/2.;
  geom->startx[3] = 0.;

  fclose(fp);

  return 0;
}

int koral_translate(double ****data, double *prims, geom_koral *geom) {
  // uses data array (koral-format defined) to populate the prims
  // array (intermediate/iharm-like format)

  // we're hard-coded to a certain coordinate configuration right now
  assert(geom->zone_coordinates == MKS3COORDS);
  assert(geom->file_coordinates == KORAL_BLCOORDS);

  #pragma omp parallel for collapse(2)
  for (int i=0; i<geom->nx; ++i) {
    fprintf(stderr, "%d ", i);
    if (i+1 == geom->nx) fprintf(stderr, "\n");
    for (int j=0; j<geom->ny; ++j) {
      for (int k=0; k<geom->nz; ++k) {

        // first get the coordinates
        double xMKS3[4] = { 0 };
        double xEKS[4] = { 0 };
        double xKS[4] = { 0 };
        double xBL[4] = { 0 };
        xMKS3[1] = geom->startx[1]+(i+0.5)*geom->dx[1];
        xMKS3[2] = geom->startx[2]+(j+0.5)*geom->dx[2];
        xMKS3[3] = geom->startx[3]+(k+0.5)*geom->dx[3];
        coco_MKS32KS(xMKS3, xKS);
        coco_KS2EKS(xKS, xEKS);
        coco_KS2BL(xKS, xBL);

        double gam = 5./3;
        double rho = data[i][j][k][3];
        double tgas = data[i][j][k][4];
        double ucon[4] = { 0. };
        double bcon[4] = { 0. };
        double te, ti;

        ucon[0] = data[i][j][k][5];
        ucon[1] = data[i][j][k][6];
        ucon[2] = data[i][j][k][7];
        ucon[3] = data[i][j][k][8];

        if (geom->dtype == KORAL_GRMHD) {

          double gcovBL[4][4] = { 0 };
          calc_gcon_bl(xBL, gcovBL);

          double ucov[4] = { 0 };
          ab_b_a(gcovBL, ucon, ucov);

          bcon[1] = data[i][j][k][10];
          bcon[2] = data[i][j][k][11];
          bcon[3] = data[i][j][k][12];

          bcon[0] = -1. * (bcon[1]*ucov[1] + bcon[2]*ucov[2] + bcon[3]*ucov[3]) / ucov[0];

          gam = geom->gam;

        } else if (geom->dtype == KORAL_RADGRMHD) {

          bcon[0] = data[i][j][k][9];
          bcon[1] = data[i][j][k][10];
          bcon[2] = data[i][j][k][11];
          bcon[3] = data[i][j][k][12];

          te = data[i][j][k][13];
          ti = data[i][j][k][14];

          gam = data[i][j][k][15];

        }

        if (0 == 1) {
          // check whether u.u == -1
          double ucovBL[4] = { 0 };
          double gconBL[4][4] = { 0 };
          double gcovBL[4][4] = { 0 };

          calc_gcon_bl(xBL, gconBL);
          invert_44(gconBL, gcovBL);

          ab_b_a(gcovBL,ucon,ucovBL);

          fprintf(stderr, "%d %d %d -> %g\n", i,j,k, a_a_(ucon,ucovBL));
        }

        // translate to prims and update prims array

        double UU = tgas * rho / ( gam - 1. );

        // the four-vectors are output in OUTCOORDS (BLCOORDS) with those
        // components so we need to translate them into the MKS3 basis

        // compute the Jacobian
        double dxdxBL2KS[4][4] = { 0 };
        double dxdxKS2EKS[4][4] = { 0 };
        dxdx_BL2KS(xBL, dxdxBL2KS);
        dxdx_KS2EKS(xKS, dxdxKS2EKS);

        // change basis
        double uconKS[4], uconEKS[4], bconKS[4], bconEKS[4];
        ab_b_a(dxdxBL2KS, ucon, uconKS);
        ab_b_a(dxdxBL2KS, bcon, bconKS);
        ab_b_a(dxdxKS2EKS, uconKS, uconEKS);
        ab_b_a(dxdxKS2EKS, bconKS, bconEKS);

        // metric &c for the next few stages
        double gconEKS[4][4] = { 0 };
        double gcovEKS[4][4] = { 0 };
        calc_gcon_eks(xEKS, gconEKS);
        invert_44(gconEKS, gcovEKS);

        // correction stage where we enforce u.u=-1 and u.b=0
        // shouldn't affect good zones
        double ucovEKS[4] = { 0 };
        ab_b_a(gcovEKS, uconEKS, ucovEKS);
        if (fabs(a_a_(uconEKS,ucovEKS)+1.) > 1.e-1) {
          set_u0(uconEKS, gconEKS, gcovEKS);
          ab_b_a(gcovEKS, uconEKS, ucovEKS);
        }
        if (fabs(a_a_(bconEKS,ucovEKS)) > 1.e-1) {
          bconEKS[0] = -1. * (bconEKS[1]*ucovEKS[1] + bconEKS[2]*ucovEKS[2] + bconEKS[3]*ucovEKS[3]) / ucovEKS[0];
        }

        if (0 == 1) {
          // check, did we transform things correctly? if this doesn't 
          // return reasonable numbers, then we know we didn't...
          double gconEKS[4][4] = { 0 };
          double gcovEKS[4][4] = { 0 };
          calc_gcon_eks(xEKS, gconEKS);
          invert_44(gconEKS, gcovEKS);

          double ucovEKS[4], bcovEKS[4];
          ab_b_a(gcovEKS,uconEKS,ucovEKS);
          ab_b_a(gcovEKS,bconEKS,bcovEKS);

          double gconBL[4][4] = { 0 };
          double gcovBL[4][4] = { 0 };
          calc_gcon_bl(xBL, gconBL);
          invert_44(gconBL, gcovBL);

          double ucovBL[4], bcovBL[4];
          ab_b_a(gcovBL,ucon,ucovBL);
          ab_b_a(gcovBL,bcon,bcovBL);

          fprintf(stderr, "%d %d %d -> (EKS) %g %g (BL) %g %g\n", 
            i,j,k, a_a_(uconEKS,ucovEKS), a_a_(uconEKS,bcovEKS),
            a_a_(ucon,ucovBL), a_a_(ucon,bcovBL));

          exit(41);
        }

        // finally get the primitives from uconEKS, bconEKS

        // first get v^i
        double alpha = 1. / sqrt(-gconEKS[0][0]);
        double gamma = uconEKS[0] * alpha;
        double v1 = uconEKS[1] + gamma * alpha * gconEKS[0][1];
        double v2 = uconEKS[2] + gamma * alpha * gconEKS[0][2];
        double v3 = uconEKS[3] + gamma * alpha * gconEKS[0][3];

        // then get B^i
        double B1 = bconEKS[1]*uconEKS[0] - bconEKS[0]*uconEKS[1];
        double B2 = bconEKS[2]*uconEKS[0] - bconEKS[0]*uconEKS[2];
        double B3 = bconEKS[3]*uconEKS[0] - bconEKS[0]*uconEKS[3];

        // Find thetae from temperature in K, and entropy from that
        double thetae = te*K_BOLTZ_CGS/(M_ELECTRON_CGS * CL_CGS * CL_CGS);
        double Thetae_unit = (M_PROTON_CGS / M_ELECTRON_CGS);
        double kel = thetae/pow(rho * geom->rhoCGS2GU , geom->gam_e-1.)/Thetae_unit;
        double ktot = 0.;  // TODO Add if this becomes necessary

        // update prims array
        int nprim;
        if (geom->dtype == KORAL_RADGRMHD) {nprim = 10;} else {nprim = 8;}
        prims[IJKP(i,j,k,0,geom->ny,geom->nz,nprim)] = rho * geom->rhoCGS2GU;
        prims[IJKP(i,j,k,1,geom->ny,geom->nz,nprim)] = UU * geom->uintCGS2GU;
        prims[IJKP(i,j,k,2,geom->ny,geom->nz,nprim)] = v1;
        prims[IJKP(i,j,k,3,geom->ny,geom->nz,nprim)] = v2;
        prims[IJKP(i,j,k,4,geom->ny,geom->nz,nprim)] = v3;
        prims[IJKP(i,j,k,5,geom->ny,geom->nz,nprim)] = B1 * sqrt(geom->endenCGS2GU);
        prims[IJKP(i,j,k,6,geom->ny,geom->nz,nprim)] = B2 * sqrt(geom->endenCGS2GU);
        prims[IJKP(i,j,k,7,geom->ny,geom->nz,nprim)] = B3 * sqrt(geom->endenCGS2GU);
        if (geom->dtype == KORAL_RADGRMHD) {
          prims[IJKP(i,j,k,8,geom->ny,geom->nz,nprim)] = kel;
          prims[IJKP(i,j,k,9,geom->ny,geom->nz,nprim)] = ktot;
        }

        if (0 == 1) {
          // did we do that correctly? compute ucon, bcon in KORAL's 
          // BLCOORDS and compare to the original values...

          double Vcon[4] = { 0., v1, v2, v3 };
          double Ucon[4] = { 0 }, Ucov[4] = { 0 };
          double Bcon[4] = { 0 };

          // start from scratch and get coordinates / metric
          double x[4] = { 0 };
          double xks[4] = { 0 };
          double xeks[4] = { 0 };
          double gcon[4][4] = { 0 }, gcov[4][4] = { 0 };
          x[1] = geom->startx[1]+(i+0.5)*geom->dx[1];
          x[2] = geom->startx[2]+(i+0.5)*geom->dx[2];
          x[3] = geom->startx[3]+(i+0.5)*geom->dx[3];
          coco_MKS32KS(x, xks);
          coco_KS2EKS(xks, xeks);
          calc_gcon_eks(xeks, gcon);
          invert_44(gcon, gcov);

          // get four-velocity
          double VdotV = 0.;
          for (int i = 1; i < 4; ++i) {
            for (int j = 1; j < 4; ++j) {
              VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
            }
          }
          double Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
          Ucon[0] = -Vfac * gcon[0][0];
          for (int i = 1; i < 4; ++i) {
            Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
          }
          ab_b_a(gcov, Ucon, Ucov);

          // get bcon
          Bcon[0] = B1*Ucov[1] + B2*Ucov[2] + B3*Ucov[3];
          Bcon[1] = (B1 + Ucon[1] * Bcon[0]) / Ucon[0];
          Bcon[2] = (B2 + Ucon[2] * Bcon[0]) / Ucon[0];
          Bcon[3] = (B3 + Ucon[3] * Bcon[0]) / Ucon[0];

          // the original "prims" were in BL coords, so get the Jacobians
          // back to there
          double dxdxEKS2KS[4][4] = { 0 };
          double dxdxKS2BL[4][4] = { 0 };
          dxdx_EKS2KS(xeks, dxdxEKS2KS);
          dxdx_KS2BL(xks, dxdxKS2BL);

          double UconKS[4] = { 0 }, UconBL[4] = { 0 };
          double BconKS[4] = { 0 }, BconBL[4] = { 0 };
          ab_b_a(dxdxEKS2KS, Ucon, UconKS);
          ab_b_a(dxdxEKS2KS, Bcon, BconKS);
          ab_b_a(dxdxKS2BL, UconKS, UconBL);
          ab_b_a(dxdxKS2BL, BconKS, BconBL);

          P4VEC("Ucon", UconBL);
          fprintf(stderr, "Ucon %g %g %g %g\n", ucon[0], ucon[1], ucon[2], ucon[3]);

          P4VEC("Bcon", BconBL);
          fprintf(stderr, "Bcon %g %g %g %g\n", bcon[0], bcon[1], bcon[2], bcon[3]);

          exit(43);
        }

      }
    }
  }

  return 0;
}

void set_u0(double ucon[4], double gcon[4][4], double gcov[4][4]) {
  // enforce u.u = -1 by setting ucon[0]
  // if we can't, then set to normal observer velocity

  double a = gcov[0][0];
  double b = 0.;
  double c = 1.;

  for (int i=1; i<4; ++i) {
    b += 2. * ucon[i] * gcov[0][i];
    for (int j=1; j<4; ++j) {
      c += ucon[i] * ucon[j] * gcov[i][j];
    }
  }

  double delta = b*b - 4.*a*c;

  if (delta < 0.) {

    double ucov[4] = { 0., 0., 0., 0. };
    ucov[0] = sqrt(-1./gcon[0][0]);
    ab_b_a(gcon, ucov, ucon);

    return;
  }

  if (a < 0.) {
    ucon[0] = (-b - sqrt(delta)) / 2. / a;
  } else {
    ucon[0] = (-b + sqrt(delta)) / 2. / a;
  }

}
