#include "utils.h"
#include "metrics.h"
#include "coordinates.h"

#include "iharm/iharm.h"
#include "koral/koral.h"

int OUTPUT_GRID = 0;
int OUTPUT_GRID_ONLY = 0;

void set_ofname(char *fname, char ofname[256]) {
  // replaces .dat with .h5
  char *s = fname, *e = ofname;
  char *lastdot = 0;
  while (*s != 0) {
    if (*s == '.') lastdot = e;
    *e++ = *s++;
  }
  *lastdot++ = '.';
  *lastdot++ = 'h';
  *lastdot++ = '5';
  *lastdot = 0;
}

// $ ./babel datafile defines
int main(int argc, char *argv[]) {

  // do a safety check to make sure we won't overwrite anything at dump time
  char ofname[256];
  set_ofname(argv[2], ofname);
  FILE *ofp = fopen(ofname, "r");
  if (ofp) {
    fprintf(stderr, " ! refusing to continue. file exists at %s\n", ofname);
    exit(-1);
  }

  geom_koral geom = { 0 };
  //geom.dtype = KORAL_GRMHD;
  geom.dtype = KORAL_RADGRMHD;

  fprintf(stderr, "reading KORAL define for geometry and sizes...\n");
  if (koral_init(argv[1], &geom) != 0) {
    fprintf(stderr, " ! unable to get parameters. exiting.\n");
    exit(-1);
  }

  // TODO allocate GRMHD/Rad differently here


  double ****rawdata = malloc4(geom.nx, geom.ny, geom.nz, geom.np);
  int nprim;
  if (geom.dtype == KORAL_RADGRMHD) {nprim = 10;} else {nprim = 8;}
  double *prims = malloc(sizeof(*prims) * geom.nx * geom.ny * geom.nz * nprim);

  fprintf(stderr, "reading KORAL datafile...\n");
  if (koral_read(argv[2], rawdata, &geom) != 0) {
    fprintf(stderr, " ! unable to read KORAL datafile. exiting.\n");
    exit(-2);
  }
  fprintf(stderr, "  n[i]      = (%d, %d, %d)\n", geom.nx, geom.ny, geom.nz);
  fprintf(stderr, "  startx[i] = (%g, %g, %g)\n", geom.startx[1], geom.startx[2], geom.startx[3]);
  fprintf(stderr, "  dx[i]     = (%g, %g, %g)\n", geom.dx[1], geom.dx[2], geom.dx[3]);

  if (OUTPUT_GRID || OUTPUT_GRID_ONLY) {
    iharm_write_grid_koral("grid.h5", &geom);
  }

  if (! OUTPUT_GRID_ONLY) {
    fprintf(stderr, "transforming data to iharm3d format...\n");
    koral_translate(rawdata, prims, &geom);

    fprintf(stderr, "saving iharm dump file...\n");
    iharm_write_header_koral(ofname, &geom);
    iharm_dump(ofname, prims, geom.nx, geom.ny, geom.nz, nprim);
  }

  fprintf(stderr, "done!\n\n");

  // housekeeping
  /*  -- the OS is smart enough to reclaim the memory. don't waste time with iterations.
  free(prims);
  free4(geom.nx, geom.ny, geom.nz, geom.np, rawdata);
   */

}
