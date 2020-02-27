#!/usr/bin/env python3

## Translator of BHAC (binary) output to HARM format
# Currently only supports output in LOG_KS and BHAC_MKS metrics, but is extensible

# Usage: bhac-translate.py file1.blk [file2.blk ...]

import sys
import numpy as np
import h5py

bhacfnames = sys.argv[1:]
harmfnames = [bfname[:-3]+"h5" for bfname in bhacfnames]

for bfname, harmfname in zip(bhacfnames, harmfnames):
  # Open bhac fluid dump
  with open(bfname,"rb") as bf:
    N1 = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    N2 = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    N3 = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    t = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    a = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    startx1 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    startx2 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    startx3 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    stopx1 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    stopx2 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    stopx3 = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    gam = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    hslope = np.frombuffer(bf.read(8), dtype=np.float64)[0]
    # metric is a flag to the ray-tracing code for the co-ordinate system in which the GRMHD data is given: -1 is Minkowski, 0 is Boyer-Lindquist, 1 is Kerr-Schild, 2 is modified Kerr-Schild (as used in HARM), 3 is logarithmic Kerr-Schild (i.e. Kerr-Schild with a logarithmic radial co-ordinate), 4 is Cartesian Kerr-Schild.
    # Note: Modern BHAC MKS only differs from HARM MKS in using X2 in [0,Pi] and hslope->(1-hslope) i.e. LOG_KS is hslope 0
    # TODO it wouldn't be too hard to just go ahead and convert everything here, writing in Gammie MKS
    bhac_metric = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    metric = [ "MINKOWSKI", "BL", "KS", "BHAC_MKS", "LOG_KS", "CART_KS" ][bhac_metric+1]
    # code is a flag to the ray-tracing code, equal to 0 for input data from BHAC, equal to 1 for input data from RAISHIN, and equal to 2 for input data from HARM.
    code = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    # dim is a flag to the ray-tracing code for the dimenisonality and form of the data: 2/3 for 2D/3D, followed by 12 for  data or 13 for  data (2D only). Consequently 2D  data means dim = 212 and 2D  data means dim = 213. All 3D data is specified by dim = 3.
    dim = np.frombuffer(bf.read(4), dtype=np.int32)[0]
    # prec is a (logical) flag to the ray-tracing code for the numerical precision of the GRMHD data (not the header). prec = .TRUE. implies single precision data, whereas prec = .FALSE. implies double precision data.
    prec = np.frombuffer(bf.read(4), dtype=np.int32)[0]

    # Because apparently headers need padding idk
    bf.seek(256, 0)

    print("Time: {} File size: {}x{}x{}, Spin: {}, hslope: {}, Grid start: ({},{},{}) stop: ({},{},{})".format(
           t, N1, N2, N3, a, hslope, startx1, startx2, startx3, stopx1, stopx2, stopx3))

    prims_bhac = np.zeros((N1,N2,N3,8), dtype=np.float32)
    for k in range(N3):
      for j in range(N2):
        for i in range(N1):
          prims_bhac[i,j,k,:] = np.frombuffer(bf.read(8*4), dtype=np.float32, count=8)

  # Now write it all
  with h5py.File(harmfname,"w") as hf:
    hf['t'] = t
    hf['dump_cadence'] = 10 # TODO is there any way to read this from BHAC output?
    # Basic header flags
    hdr = hf.create_group('header')
    hdr['n1'] = N1
    hdr['n2'] = N2
    hdr['n3'] = N3
    hdr['gam'] = gam
    hdr['has_electrons'] = False
    hdr['n_prim'] = 8
    hdr['n_prims_passive'] = 0
    # We'll output in the HARM order/values
    prim_names = [ b"RHO", b"UU", b"U1", b"U2", b"U3", b"B1", b"B2", b"B3" ]
    hdr.create_dataset("prim_names", data=np.array(prim_names, dtype='S'))

    # Non-coordinate-specific geometry
    geom = hf.create_group('header/geom')
    geom['startx1'] = startx1
    geom['startx2'] = startx2
    geom['startx3'] = startx3
    geom['dx1'] = (stopx1 - startx1)/N1
    geom['dx2'] = (stopx2 - startx2)/N2
    geom['dx3'] = (stopx3 - startx3)/N3
    geom['n_dim'] = dim + 1 # HARM includes _all_ dimensions

    hf.create_dataset("/header/metric", data=np.string_("BHAC_MKS"))
    bmks = hf.create_group('header/geom/bhac_mks')
    # Patch a typo
    if a > 0.98:
      bmks['a'] = 0.96875
    else:
      bmks['a'] = a

    bmks['r_in'] = np.exp(startx1)
    bmks['r_out'] = np.exp(stopx1)

    if metric == "BHAC_MKS":
      bmks['hslope'] = hslope
    elif metric == "LOG_KS":
      bmks['hslope'] = 0.
    else:
      print("Metric {} not supported!".format(metric))
      exit(-1)

    # Convert primitives
    prims = np.zeros_like(prims_bhac)
    # RHO is the same
    prims[:,:,:,0] = prims_bhac[:,:,:,0]
    # UU from P
    prims[:,:,:,1] = prims_bhac[:,:,:,4]/(gam-1)
    # v^1,2,3 * gamma == u^1,2,3
    prims[:,:,:,2] = prims_bhac[:,:,:,1]
    prims[:,:,:,3] = prims_bhac[:,:,:,2]
    prims[:,:,:,4] = prims_bhac[:,:,:,3]
    # B^i defined the same way
    prims[:,:,:,5] = prims_bhac[:,:,:,5]
    prims[:,:,:,6] = prims_bhac[:,:,:,6]
    prims[:,:,:,7] = prims_bhac[:,:,:,7]

    hf['prims'] = prims
