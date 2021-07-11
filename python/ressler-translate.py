# python convert_cks.py --is=[starti] --ie=[endi] --file_prefix=star_wind.out2
# python translate_npz_iharm.py [athinput.gr_init] [npz-dump-file]

import numpy as np
import time
import h5py
import sys
import os

def get_gcov_ks(R, H, a=0):
    # calculate for ks components in ks coordinates at r,h
    gcov = np.zeros((R.shape[0],R.shape[1],4,4))
    cth = np.cos(H)
    sth = np.sin(H)
    s2 = sth*sth
    rho2 = R*R + a*a*cth*cth
    gcov[:,:,0,0] = (-1. + 2. * R / rho2)
    gcov[:,:,0,1] = (2. * R / rho2)
    gcov[:,:,0,3] = (-2. * a * R * s2 / rho2)
    gcov[:,:,1,0] = gcov[:,:,0,1]
    gcov[:,:,1,1] =  (1. + 2. * R / rho2)
    gcov[:,:,1,3] =  (-a * s2 * (1. + 2. * R / rho2))
    gcov[:,:,2,2] = rho2
    gcov[:,:,3,0] = gcov[:,:,0,3]
    gcov[:,:,3,1] = gcov[:,:,1,3]
    gcov[:,:,3,3] = s2 * (rho2 + a*a * s2 * (1. + 2. * R / rho2))
    return gcov

def get_gcon_eks_3d(R, H, a=0):
    gcov_ks = get_gcov_ks(R[:, :, 0], H[:, :, 0], a=a)
    gcon_ks = np.linalg.inv(gcov_ks)
    dXdx = np.zeros((N1, N2, 4, 4))  # X is eks, x is ks
    dXdx[:, :, 0, 0] = 1.
    dXdx[:, :, 1, 1] = 1./R[:, :, 0]
    dXdx[:, :, 2, 2] = 1.
    dXdx[:, :, 3, 3] = 1.
    gcon_eks = np.einsum('abki,abkj->abij',dXdx,np.einsum('ablj,abkl->abkj',dXdx,gcon_ks))
    gcon_eks_3d = np.zeros((*R.shape, 4, 4))
    gcon_eks_3d[:,:,:,:,:] = gcon_eks[:,:,None,:,:]
    return gcon_eks_3d

def write_header(hfp, r, h, p, gam, a, mu_i, mu_e, mu_tot):
    N1, N2, N3 = r.shape
    
    hfp.create_group('header')
    hfp.create_group('header/weights')
    hfp.create_group('header/geom')
    hfp.create_group('header/geom/eks')

    hfp.create_dataset("/header/metric", data=np.string_("EKS"))
    hfp['header']['n1'] = N1
    hfp['header']['n2'] = N2
    hfp['header']['n3'] = N3
    
    dx1 = np.diff(np.log(r[:, 0, 0])).mean()
    dx2 = np.diff(h[0, :, 0]).mean()
    dx3 = np.diff(p[0, 0, :]).mean()
    startx1 = np.log(r[0, 0, 0]) - dx1/2.
    startx2 = h[0, 0, 0] - dx2/2.
    startx3 = 0.

    # output weights
    hfp['header']['weights']['mu_i'] = mu_i
    hfp['header']['weights']['mu_e'] = mu_e
    hfp['header']['weights']['mu_tot'] = mu_tot

    # this is the left edge of the grid
    hfp['header']['geom']['startx1'] = startx1
    hfp['header']['geom']['startx2'] = startx2
    hfp['header']['geom']['startx3'] = startx3
    
    # this is the separation between grid zone centers
    hfp['header']['geom']['dx1'] = dx1
    hfp['header']['geom']['dx2'] = dx2
    hfp['header']['geom']['dx3'] = dx3

    hfp['header']['geom']['eks']['a'] = a

    # these give the actual boundaries
    hfp['header']['geom']['eks']['r_eh'] = 1. + np.sqrt(1. - a*a)
    hfp['header']['geom']['eks']['r_in'] = np.exp(startx1)
    hfp['header']['geom']['eks']['r_out'] = np.exp(startx1 + dx1*N1)
    
    hfp['header']['n_prim'] = 2+3+3+1
    hfp['header']['n_prims_passive'] = 0
    hfp['header']['gam'] = gam
    hfp['header']['has_electrons'] = 3  # forces Theate_unit = MP/ME
    hfp['header']['has_radiation'] = 0
    
    prim_names = [ b"RHO", b"UU", b"U1", b"U2", b"U3", b"B1", b"B2", b"B3" ]
    hfp['header'].create_dataset("prim_names", data=np.array(prim_names, dtype='S'))


def get_from_athinput(fname, skey):
    fp = open(fname, 'r')
    for line in fp.readlines():
        if "=" not in line:
            continue
        try:
            line = line.replace("=", "").split()
            key = line[0]
            val = line[1]
            if key == skey:
                return val
        except:
            pass


if __name__ == "__main__":

    # deal with units here. by convention, we choose that
    # the M_unit input to ipole/grmonty should be one, so
    # we need to scale Ne, B, and Thetae appropriately. I
    # think the easy way to do this is to recover the cgs
    # units here first, then back out the values required
    # to work with M_unit = 1. we assume that we consider
    # SgrA*.

    # THESE WILL BE SET BY IPOLE:
    # HERE WE ASSUME MBH = 4.14e6
    # Ne_unit = (M_unit/L_unit^3) / (MP + ME)
    # B_unit = CL * Sqrt[ 4 Pi M_unit / L_unit^3 ]

    if True:
        # beta_w = 1e2
        model = {
          'rho_unit': 53939.4,
          'P_unit': 5.05615e+09,
          'B_unit': 71106.6
        }

    if False:
        # beta_w = 1e6
        model = {
            'rho_unit': 29737.4,
            'P_unit': 2.78751e+09,
            'B_unit': 52796.9
        }

    # compute target units in ipole
    pc = 3.086e18
    kyr = 3.154e10
    MBH = 4.3e6
    GNEWT = 6.6742e-8
    CL = 2.99792458e10
    MSUN = 1.989e33
    MP = 1.67262171e-24
    ME = 9.1093826e-28
    KBOL = 1.3806505e-16
    GRAMS_TO_MeV = 5.609588357187173e+26
    keV_to_Kelvin = 1.16045e7
    mp_over_kev = 9.994827

    # these are the units ipole will use
    target_L_unit = MBH * MSUN * GNEWT / (CL * CL)
    target_Ne_unit = (1. / target_L_unit**3) / (MP + ME)
    target_B_unit = CL * np.sqrt(4. * np.pi / target_L_unit**3)

    # weights
    mu_i = 4.16
    mu_e = 2.0
    mu_tot = 1.351

    # target rho conversion factor
    rho_to_n_cgs = 40.46336
    GRMHD_RHO_TO_NE = rho_to_n_cgs/mu_e * model['rho_unit']
    RHO_FACTOR = GRMHD_RHO_TO_NE / target_Ne_unit * mu_e

    # target B conversion factor
    GRMHD_B_TO_B = model['B_unit'] * pc/kyr * np.sqrt(4.*np.pi * MSUN/pc**3)
    B_FACTOR = GRMHD_B_TO_B / target_B_unit

    # target Thetae conversion factor
    TE_FACTOR = model['P_unit'] / model['rho_unit'] * mu_tot * mp_over_kev * keV_to_Kelvin
    KELVIN_TO_THETAE = KBOL / ME / CL / CL

    # begin conversion
    t0 = time.time()

    overwrite = "--overwrite" in sys.argv

    fnames = [x for x in sys.argv if x[-4:]==".npz"]
    if len(fnames) != 1:
      print("! must pass one .npz input file. halting.")
      exit(1)
    fname = fnames[0]

    ofname = fname.replace("npz", "h5")
    if os.path.exists(ofname) and not overwrite:
      print("! output file {0:s} exists. rerun with --overwrite?".format(ofname))
      exit(2)

    # load from parameter file
    athinput = sys.argv[1]
    fluid_gamma = float(get_from_athinput(athinput, "gamma"))
    bhspin = float(get_from_athinput(athinput, "a"))
    dump_cadence = float(get_from_athinput(athinput, "dt"))  ## WARNING: assumes all dt are the same

    # load data from file
    print(f" - loading {fname}")
    d = np.load(fname)
    t = d['t']
    r_ks = d['r']
    h_ks = d['th']
    p_ks = d['ph']
    x = d['x']
    y = d['y']
    z = d['z']
    bsq = d['bsq']
    rho = d['rho']
    press = d['press']
    ucon_ks = d['uu'].transpose((1,2,3,0))
    bcon_ks = d['bu'].transpose((1,2,3,0))
    N1, N2, N3 = rho.shape
    print(" - got {0:d}x{1:d}x{2:d} file at t={3:g}".format(N1, N2, N3, t))

    # perform validation
    print(" - validating u.u and b.b in 2d")
    gcov_ks = get_gcov_ks(r_ks[:, :, 0], h_ks[:, :, 0], a=bhspin)
    ucov_ks = np.einsum('abij,abi->abj', gcov_ks, ucon_ks[:, :, 0, :])
    bcov_ks = np.einsum('abij,abi->abj', gcov_ks, bcon_ks[:, :, 0, :])
    udotu = np.einsum('abi,abi->ab', ucon_ks[:, :, 0, :], ucov_ks)
    bsq_check = np.einsum('abi,abi->ab', bcon_ks[:, :, 0, :], bcov_ks)
    print("   udotu check:", np.allclose(udotu, -1.))
    print("   bsq check:", np.allclose(bsq_check, bsq[:, :, 0]))

    # note here that X is eks and x is ks
    print(" - translating from native ks -> eks")
    dXdx = np.zeros((N1, N2, N3, 4, 4))
    dXdx[:, :, :, 0, 0] = 1.
    dXdx[:, :, :, 1, 1] = 1./r_ks
    dXdx[:, :, :, 2, 2] = 1.
    dXdx[:, :, :, 3, 3] = 1.
    ucon_eks = np.einsum('abciu,abci->abcu', dXdx, ucon_ks)
    bcon_eks = np.einsum('abciu,abci->abcu', dXdx, bcon_ks)
   
    # generate prims with eks components at eks spacing
    print(" - generating prims array")
    prims = np.zeros((N1, N2, N3, 9))
    prims[:, :, :, 0] = rho * RHO_FACTOR
    prims[:, :, :, 1] = press / (fluid_gamma - 1.) * B_FACTOR*B_FACTOR

    gcon_eks = get_gcon_eks_3d(r_ks, h_ks, a=bhspin)

    prims[:, :, :, 2] = ucon_eks[:, :, :, 1] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 1] / gcon_eks[:, :, :, 0, 0]
    prims[:, :, :, 3] = ucon_eks[:, :, :, 2] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 2] / gcon_eks[:, :, :, 0, 0]
    prims[:, :, :, 4] = ucon_eks[:, :, :, 3] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 3] / gcon_eks[:, :, :, 0, 0]
    
    prims[:, :, :, 5] = bcon_eks[:, :, :, 1] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 1]
    prims[:, :, :, 6] = bcon_eks[:, :, :, 2] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 2]
    prims[:, :, :, 7] = bcon_eks[:, :, :, 3] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 3]

    prims[:, :, :, 5:8] *= B_FACTOR

    prims[:, :, :, 8] = press / rho * TE_FACTOR * KELVIN_TO_THETAE

    # output
    print(" - writing to {0:s}".format(ofname))
    hfp = h5py.File(ofname, 'w')
    write_header(hfp, r_ks, h_ks, p_ks, fluid_gamma, bhspin, mu_i, mu_e, mu_tot)
    hfp['prims'] = prims
    hfp['t'] = t
    hfp['dump_cadence'] = dump_cadence
    hfp.close()

    print(" - {0:.1f} s elapsed".format(time.time() - t0))

