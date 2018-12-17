import h5py
import numpy as np
import sys

if __name__ == "__main__":

 fnames = sys.argv[1:]

 for fname in fnames:

   print("translating {0:s}".format(fname))
   ofname = fname + ".h5"

   # parse header
   fp = open(fname,'r')
   header = fp.readline().split()
   fp.close()

   t = float(header[0])
   N1 = int(header[1])
   N2 = int(header[2])
   N3 = int(header[3])
   startx1 = float(header[4])
   startx2 = float(header[5])
   startx3 = float(header[6])
   dx1 = float(header[7])
   dx2 = float(header[8])
   dx3 = float(header[9])

   MBH = float(header[12])
   a = float(header[13])

   Lunit = float(header[14])
   Tunit = float(header[15])
   Munit = float(header[16])
   print(Lunit, Tunit, Munit)

   gam = float(header[18])
   game = float(header[19])
   gamp = float(header[20])

   rin = float(header[34])
   rout = float(header[35])
   hslope = float(header[36])
   R0 = float(header[37])

   # get prims data
   data = np.loadtxt(fname, skiprows=1)

   # get output file
   hfp = h5py.File(ofname, 'w')

   ## write header
   hfp.create_group("/header")
   hfp.create_group("/header/geom")
   hfp.create_group("/header/geom/mks")

   hfp['header']['has_electrons'] = 1

   hfp['header']['n1'] = N1
   hfp['header']['n2'] = N2
   hfp['header']['n3'] = N3
   hfp['header']['n_prim'] = 10

   hfp['header']['geom']['startx1'] = startx1
   hfp['header']['geom']['startx2'] = startx2
   hfp['header']['geom']['startx3'] = startx3
   hfp['header']['geom']['dx1'] = dx1
   hfp['header']['geom']['dx2'] = dx2
   hfp['header']['geom']['dx3'] = dx3

   hfp.create_dataset("/header/metric", data=np.string_("MKS"))
   hfp['header']['geom']['mks']['a'] = a
   hfp['header']['geom']['mks']['Rin'] = rin
   hfp['header']['geom']['mks']['Rout'] = rout
   hfp['header']['geom']['mks']['R0'] = R0
   hfp['header']['geom']['mks']['hslope'] = hslope

   hfp['header']['gam'] = gam
   hfp['header']['gam_e'] = game
   hfp['header']['gam_p'] = gamp

   hfp['t'] = t

   ## write prims
   prims = np.zeros((N1, N2, N3, 10))
   nprim = 0
   for vals in data:
     j = nprim % N2
     i = (nprim - j) // N2
     k = 0
     nprim += 1
     prims[i,j,k,0] = vals[4]
     prims[i,j,k,1] = vals[5]
     prims[i,j,k,2] = vals[6]
     prims[i,j,k,3] = vals[7]
     prims[i,j,k,4] = vals[8]
     prims[i,j,k,5] = vals[9]
     prims[i,j,k,6] = vals[10]
     prims[i,j,k,7] = vals[11]
     prims[i,j,k,8] = vals[13]
     prims[i,j,k,9] = vals[12]

   hfp['prims'] = prims
   hfp.close()
