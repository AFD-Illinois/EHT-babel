
#include "iharm.h"

// write header from koral geometry
int iharm_write_header_koral(char *fname, geom_koral *geom) {

  assert(geom->zone_coordinates == MKS3COORDS); 

  hid_t fid = H5Fcreate(fname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) {
    return -1;
  }

  // create groups
  h5io_add_group(fid, "/header");
  h5io_add_group(fid, "/header/geom");
  h5io_add_group(fid, "/header/geom/mks3"); // TODO don't hardcode this

  h5io_add_data_int(fid, "/header/n1", geom->nx);
  h5io_add_data_int(fid, "/header/n2", geom->ny);
  h5io_add_data_int(fid, "/header/n3", geom->nz);
  int nprim, has_elec;
  if (geom->dtype == KORAL_RADGRMHD) {
    nprim = 10;
    has_elec = 1;
  } else {
    nprim = 8;
    has_elec = 0;
  }
  h5io_add_data_int(fid, "/header/n_prim", nprim);
  h5io_add_data_int(fid, "/header/has_electrons", has_elec);

  h5io_add_data_dbl(fid, "/header/geom/startx1", geom->startx[1]);
  h5io_add_data_dbl(fid, "/header/geom/startx2", geom->startx[2]);
  h5io_add_data_dbl(fid, "/header/geom/startx3", geom->startx[3]);
  h5io_add_data_dbl(fid, "/header/geom/dx1", geom->dx[1]);
  h5io_add_data_dbl(fid, "/header/geom/dx2", geom->dx[2]);
  h5io_add_data_dbl(fid, "/header/geom/dx3", geom->dx[3]);

  h5io_add_data_str(fid, "/header/metric", "MKS3"); // TODO don't hardcode?
  h5io_add_data_dbl(fid, "/header/geom/mks3/a", geom->a);
  h5io_add_data_dbl(fid, "/header/geom/mks3/R0", geom->R0);
  h5io_add_data_dbl(fid, "/header/geom/mks3/H0", geom->H0);
  h5io_add_data_dbl(fid, "/header/geom/mks3/MY1", geom->MY1);
  h5io_add_data_dbl(fid, "/header/geom/mks3/MY2", geom->MY2);
  h5io_add_data_dbl(fid, "/header/geom/mks3/MP0", geom->MP0);

  h5io_add_data_dbl(fid, "/header/gam", geom->gam);
  h5io_add_data_dbl(fid, "/header/gam_e", geom->gam_e);
  h5io_add_data_dbl(fid, "/header/gam_p", geom->gam_p);
  
  h5io_add_data_dbl(fid, "/t", geom->t);
  h5io_add_data_dbl(fid, "/dump_cadence", geom->DTOUT1);

  H5Fclose(fid);

  return 0;
}

// write harm grid file from koral geometry
int iharm_write_grid_koral(char *fname, geom_koral *geom) {

  // file coordinates
  double *X1 = malloc(sizeof(*X1) * geom->nx * geom->ny * geom->nz * 1);
  double *X2 = malloc(sizeof(*X2) * geom->nx * geom->ny * geom->nz * 1);
  double *X3 = malloc(sizeof(*X3) * geom->nx * geom->ny * geom->nz * 1);

  // KS coordinates
  double *r = malloc(sizeof(*r) * geom->nx * geom->ny * geom->nz * 1);
  double *th = malloc(sizeof(*th) * geom->nx * geom->ny * geom->nz * 1);
  double *phi = malloc(sizeof(*phi) * geom->nx * geom->ny * geom->nz * 1);

  // plotting coordinates
  double *X = malloc(sizeof(*X) * geom->nx * geom->ny * geom->nz * 1);
  double *Y = malloc(sizeof(*Y) * geom->nx * geom->ny * geom->nz * 1);
  double *Z = malloc(sizeof(*Z) * geom->nx * geom->ny * geom->nz * 1);

  // geometry
  double *gcon = malloc(sizeof(*gcon) * geom->nx * geom->ny * geom->nz * 16);
  double *gcov = malloc(sizeof(*gcov) * geom->nx * geom->ny * geom->nz * 16);
  double *gdet = malloc(sizeof(*gdet) * geom->nx * geom->ny * geom->nz * 1);
  double *lapse = malloc(sizeof(*lapse) * geom->nx * geom->ny * geom->nz * 1);

  // compute
  #pragma omp parallel for collapse(2)
  for (int i=0; i<geom->nx; ++i) {
    fprintf(stderr, "%d ", i);
    if (i+1 == geom->nx) fprintf(stderr, "\n");
    for (int j=0; j<geom->ny; ++j) {
      for (int k=0; k<geom->nz; ++k) {

        int index3d = IJKP(i,j,k,0,geom->ny,geom->nz,1);
        int index5d = IJKP(i,j,k,0,geom->ny,geom->nz,16);

        // first get and set the coordinates
        double xMKS3[4] = { 0 };
        double xEKS[4] = { 0 };
        double xKS[4] = { 0 };
        double xCART[3] = { 0 };
        xMKS3[1] = geom->startx[1]+(i+0.5)*geom->dx[1];
        xMKS3[2] = geom->startx[2]+(j+0.5)*geom->dx[2];
        xMKS3[3] = geom->startx[3]+(k+0.5)*geom->dx[3];
        coco_MKS32KS(xMKS3, xKS);
        coco_KS2EKS(xKS, xEKS);
        coco_SPH2CART(xKS+1, xCART);

        X1[index3d] = xEKS[1];
        X2[index3d] = xEKS[2];
        X3[index3d] = xEKS[3];
        r[index3d] = xKS[1];
        th[index3d] = xKS[2];
        phi[index3d] = xKS[3];
        X[index3d] = xCART[0];
        Y[index3d] = xCART[1];
        Z[index3d] = xCART[2];

        // now get various geometric values
        double gconEKS[4][4] = { 0 };
        double gcovEKS[4][4] = { 0 };
        calc_gcon_eks(xEKS, gconEKS);
        invert_44(gconEKS, gcovEKS);

        for (int mu=0; mu<4; ++mu) {
          for (int nu=0; nu<4; ++nu) {
            int offset = mu*4 + nu;
            gcon[index5d+offset] = gconEKS[mu][nu];
            gcov[index5d+offset] = gcovEKS[mu][nu];
          }
        }

        lapse[index3d] = 1./sqrt(-gconEKS[0][0]);
        gdet[index3d] = gdet_func(gcovEKS);

      }
    }
  }

  // now save  
  hid_t fid = H5Fcreate(fname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) {
    return -1;
  }

  h5io_add_data_dbl_3d(fid, "/X1", geom->nx, geom->ny, geom->nz, X1);
  h5io_add_data_dbl_3d(fid, "/X2", geom->nx, geom->ny, geom->nz, X2);
  h5io_add_data_dbl_3d(fid, "/X3", geom->nx, geom->ny, geom->nz, X3);

  h5io_add_data_dbl_3d(fid, "/r", geom->nx, geom->ny, geom->nz, r);
  h5io_add_data_dbl_3d(fid, "/th", geom->nx, geom->ny, geom->nz, th);
  h5io_add_data_dbl_3d(fid, "/phi", geom->nx, geom->ny, geom->nz, phi);

  h5io_add_data_dbl_3d(fid, "/X", geom->nx, geom->ny, geom->nz, X);
  h5io_add_data_dbl_3d(fid, "/Y", geom->nx, geom->ny, geom->nz, Y);
  h5io_add_data_dbl_3d(fid, "/Z", geom->nx, geom->ny, geom->nz, Z);

  h5io_add_data_flt_5d(fid, "/gcon", geom->nx, geom->ny, geom->nz, 4, 4, gcon);
  h5io_add_data_flt_5d(fid, "/gcov", geom->nx, geom->ny, geom->nz, 4, 4, gcov);

  h5io_add_data_dbl_3d(fid, "/lapse", geom->nx, geom->ny, geom->nz, lapse);
  h5io_add_data_dbl_3d(fid, "/gdet", geom->nx, geom->ny, geom->nz, gdet);

  H5Fclose(fid);

  return 0;
}

// dump to file
int iharm_dump(char *fname, double *prims, int nx, int ny, int nz, int np) {

  hid_t fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
  if (fid < 0) {
    return -1;
  }

  // TODO maybe use something like this to get rid of the warning?
  //double *data[nx][ny][nz] = (int (*)[3])&prims[0];
  h5io_add_data_flt_4d(fid, "/prims", nx, ny, nz, np, prims);

  H5Fclose(fid);

  return 0;
}

