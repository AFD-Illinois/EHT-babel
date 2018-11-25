
#include "iharm.h"

// various header writers
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
  h5io_add_data_int(fid, "/header/n_prim", 8);

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

