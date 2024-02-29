"""

  Mildly modified version of the bhlight translator. Takes native output from
  iharm2d_v3 and converts it into standard hdf5 format that ipole expects.

  GNW 2024.02.29

  Usage:

  $ python iharm2d-translate.py dump0000.dat dump0001.dat ...

  to create dump0000.dat.h5, dump0001.dat.h5, etc.

  Notes:

  - Since the spin of the black hole is not stored in the dump file, it must be
    specified manually in the script below. The default value is 0.

  - The script assumes that the hslope parameter = 0.3. This is probably a good
    assumption for most simulations, but it should be kept in mind.

  - The script assumes that R0 = 0. This is probably a very good assumption for
    most simulations, but it too should be kept in mind.

"""

import h5py
import numpy as np
import sys

bhspin = 0.
hslope = 0.3
R0 = 0.

if __name__ == "__main__":

    N3 = 1

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
      startx1 = float(header[3])
      startx2 = float(header[4])
      dx1 = float(header[5])
      dx2 = float(header[6])
      tf = float(header[7])
      nstep = int(header[8])
      gam = float(header[9])
      cour = float(header[10])
      DTd = float(header[11])
      DTl = float(header[12])
      DTi = float(header[13])
      DTr = int(header[14])
      dump_cnt = int(header[15])
      image_cnt = int(header[16])
      rdump_cnt = int(header[17])
      dt = float(header[18])

      print(f' - using black hole spin = {bhspin}')
      print(f' - assuming hslope = {hslope}')
      print(f' - assuming R0 = {R0}')
      print(' - read geometry:')
      print(f'   {N1} x {N2} grid')
      print(f'   startx1: {startx1}')
      print(f'   dx1: {dx1}')

      rin = 0.98 * (1. + np.sqrt(1. - bhspin**2))
      rin = R0 + np.exp(startx1)
      rout = np.exp(N1 * dx1) * (rin - R0) + R0

      print(' - calculated geometry: ')
      print(f'   rin: {rin}')
      print(f'   rout: {rout}')
      print(f'   R0: {R0}')
    
      # get prims data
      data = np.loadtxt(fname, skiprows=1).T
      # data saved are:
      # X1, X2,
      # prims[8] (= rho, UU, u1, u2, u3, B1, B2, B3),
      # divB, ucon[4], ucov[4], bcon[4], bcov[4],
      # mhd characteristics vmin, vmax, vmin, vmax,
      # gdet, jcon[4], jcov[4]

      x1vs = data[0]  # not used
      x2vs = data[1]  # not used
      prims = data[2:10] 
      gdet = data[31]  # not used

      # get output file
      hfp = h5py.File(ofname, 'w')

      ## write header
      hfp.create_group("/header")
      hfp.create_group("/header/geom")
      hfp.create_group("/header/geom/mks")
      hfp.create_group("/header/units")

      hfp['dump_cadence'] = DTd

      hfp['header']['has_electrons'] = 0

      hfp['header']['n1'] = N1
      hfp['header']['n2'] = N2
      hfp['header']['n3'] = N3
      hfp['header']['n_prim'] = 8

      hfp['header']['geom']['startx1'] = startx1
      hfp['header']['geom']['startx2'] = startx2
      hfp['header']['geom']['startx3'] = 0.
      hfp['header']['geom']['dx1'] = dx1
      hfp['header']['geom']['dx2'] = dx2
      hfp['header']['geom']['dx3'] = 2.*np.pi

      hfp.create_dataset("/header/metric", data=np.string_("MKS"))
      hfp['header']['geom']['mks']['a'] = bhspin
      hfp['header']['geom']['mks']['r_in'] = rin
      hfp['header']['geom']['mks']['r_out'] = rout
      hfp['header']['geom']['mks']['R0'] = R0
      hfp['header']['geom']['mks']['hslope'] = hslope

      hfp['header']['gam'] = gam
      hfp['header']['gam_e'] = 4./3
      hfp['header']['gam_p'] = 5./3

      hfp['t'] = t

      ## write prims
      prims = np.zeros((N1, N2, N3, 8))
      nprim = 0
      for vals in data.T:
          j = nprim % N2
          i = (nprim - j) // N2
          k = 0
          nprim += 1
          prims[i,j,k,0] = vals[2]
          prims[i,j,k,1] = vals[3]
          prims[i,j,k,2] = vals[4]
          prims[i,j,k,3] = vals[5]
          prims[i,j,k,4] = vals[6]
          prims[i,j,k,5] = vals[7]
          prims[i,j,k,6] = vals[8]
          prims[i,j,k,7] = vals[9]
  
      hfp['prims'] = prims
      hfp.close()
