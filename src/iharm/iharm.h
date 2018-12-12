#ifndef IHARM_H
#define IHARM_H

#include <assert.h>

#include "../h5io.h"
#include "../metrics.h"
#include "../koral/koral.h" // for geom_koral

// various header writers
int iharm_write_header_koral(char *fname, geom_koral *geom);
int iharm_write_grid_koral(char *fname, geom_koral *geom);

// dump to file
int iharm_dump(char *fname, double *prims, int nx, int ny, int nz, int np);

#endif // IHARM_H