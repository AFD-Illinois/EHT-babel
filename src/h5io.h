#ifndef H5IO_H
#define H5IO_H

#define H5IO_FMT_INT (H5T_STD_I32LE)
#define H5IO_FMT_DBL (H5T_IEEE_F64LE)
#define H5IO_FMT_FLT (H5T_IEEE_F32LE)
#define H5IO_FMT_STR (H5T_C_S1)

#include <string.h>

#include <hdf5.h>

// groups
void h5io_add_group(hid_t fid, const char *path);

// blobs
typedef hid_t hdf5_blob;
void h5io_add_blob(hid_t fid, const char *path, hdf5_blob blob);

// attributes
void h5io_add_attribute_int(hid_t fid, const char *path, const char *name, int attribute);
void h5io_add_attribute_dbl(hid_t fid, const char *path, const char *name, double attribute);
void h5io_add_attribute_str(hid_t fid, const char *path, const char *name, const char *attribute);

// data
void h5io_add_data_int(hid_t fid, const char *path, int data);
void h5io_add_data_dbl(hid_t fid, const char *path, double data);
void h5io_add_data_str(hid_t fid, const char *path, const char *data);

void h5io_add_data_dbl_1d(hid_t fid, const char *path, hsize_t n1, double data[n1]);
void h5io_add_data_dbl_2d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, double data[n1][n2]);
void h5io_add_data_dbl_3d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, hsize_t n3, double data[n1][n2][n3]);

void h5io_add_data_flt_4d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, hsize_t n3, hsize_t n4, double data[n1][n2][n3][n4]);

void h5io_add_data_flt_5d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, hsize_t n3, hsize_t n4, hsize_t n5, double data[n1][n2][n3][n4][n5]);

#endif // H5IO_H