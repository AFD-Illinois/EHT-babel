import numpy as np
import h5py
import matplotlib.pyplot as plt

import sys

from scipy.interpolate import RegularGridInterpolator

import metrics

def remap_data(coords, all_coords, data_to_fill):
    x1v, x2v, x3v = coords
    all_x1s, all_x2s, all_x3s = all_coords
    mb_size = 128
    nx = ny = nz = mb_size * 3
    data = np.zeros((nx, ny, nz))
    for mbi in range(len(data_to_fill)):
        x1i = x1v[mbi].min()
        x2i = x2v[mbi].min()
        x3i = x3v[mbi].min()
        xi = np.argmax(all_x1s==x1i)
        yi = np.argmax(all_x2s==x2i)
        zi = np.argmax(all_x3s==x3i)
        xf = xi + mb_size
        yf = yi + mb_size
        zf = zi + mb_size
        data[xi:xf, yi:yf, zi:zf] = data_to_fill[mbi].transpose((2, 1, 0))
    return data

def write_header(hfp, r, h, p, gam, a):

    N1, N2, N3 = r.shape

    hfp['dump_cadence'] = 5
    hfp['t'] = 1000
    
    hfp.create_group('header')
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
    
    hfp['header']['n_prim'] = 2+3+3
    hfp['header']['n_prims_passive'] = 0
    hfp['header']['gam'] = gam
    hfp['header']['has_electrons'] = 0  # forces Theate_unit = MP/ME   ## TODO
    hfp['header']['has_radiation'] = 0
    
    prim_names = [ b"RHO", b"UU", b"U1", b"U2", b"U3", b"B1", b"B2", b"B3" ]
    hfp['header'].create_dataset("prim_names", data=np.array(prim_names, dtype='S'))

if __name__ == "__main__":

    ## input
    fname = sys.argv[1]
    interp_method = 'linear'
    bhspin = 0.9375
    fluid_gamma = 4./3
    target_n1 = 384
    target_n2 = 192
    target_n3 = 192

    ### automatic below ###
    reh = 1. + np.sqrt(1. - bhspin*bhspin)

    ## load file
    print(f" - loading {fname}")

    hfp = h5py.File(fname, 'r')
    x1v = np.array(hfp['x1v'])
    x2v = np.array(hfp['x2v'])
    x3v = np.array(hfp['x3v'])
    uov = np.array(hfp['uov'])
    LogicalLocations = np.array(hfp['LogicalLocations'])
    Levels = np.array(hfp['Levels'])
    hfp.close()

    ## get positions
    def get_all_vals_unique(xs):
        all_xs = []
        for x in xs:
            all_xs += list(x)
        return np.array(sorted(list(set(all_xs))))
    all_x1s = get_all_vals_unique(x1v)
    all_x2s = get_all_vals_unique(x2v)
    all_x3s = get_all_vals_unique(x3v)

    ## reformat data
    print(" - reformatting data")

    coords = x1v, x2v, x3v
    all_coords = all_x1s, all_x2s, all_x3s

    n1,n2,n3 = 384, 384, 384

    ucon_cks = np.zeros((4, n1, n2, n3))
    ucon_cks[0, :, :, :] = remap_data(coords, all_coords, uov[2, :, :, :]).transpose((2, 1, 0))
    ucon_cks[1, :, :, :] = remap_data(coords, all_coords, uov[3, :, :, :]).transpose((2, 1, 0))
    ucon_cks[2, :, :, :] = remap_data(coords, all_coords, uov[4, :, :, :]).transpose((2, 1, 0))
    ucon_cks[3, :, :, :] = remap_data(coords, all_coords, uov[5, :, :, :]).transpose((2, 1, 0))

    bcon_cks = np.zeros((4, n1, n2, n3))
    bcon_cks[0, :, :, :] = remap_data(coords, all_coords, uov[6, :, :, :]).transpose((2, 1, 0))
    bcon_cks[1, :, :, :] = remap_data(coords, all_coords, uov[7, :, :, :]).transpose((2, 1, 0))
    bcon_cks[2, :, :, :] = remap_data(coords, all_coords, uov[8, :, :, :]).transpose((2, 1, 0))
    bcon_cks[3, :, :, :] = remap_data(coords, all_coords, uov[9, :, :, :]).transpose((2, 1, 0))

    rho = remap_data(coords, all_coords, uov[0, :, :, :]).transpose((2, 1, 0))
    pressure = remap_data(coords, all_coords, uov[1, :, :, :]).transpose((2, 1, 0))
    internal_energy = pressure / (fluid_gamma - 1.)
    pressure_mag = remap_data(coords, all_coords, uov[10, :, :, :]).transpose((2, 1, 0))
    beta = remap_data(coords, all_coords, uov[11, :, :, :]).transpose((2, 1, 0))
    sigma = remap_data(coords, all_coords, uov[12, :, :, :]).transpose((2, 1, 0))

    ## translate four-vector ucon, bcon -> corresponding primitives
    AX1, AX2, AX3 = np.meshgrid(all_x1s, all_x2s, all_x3s, indexing='ij')
    gcon_cks = metrics.cks_inverse_metric(AX1, AX2, AX3, bhspin)

    print(" - reconstructing U primitives")
    Uprims = np.zeros((n1, n2, n3, 3))
    Uprims[:, :, :, 0] = ucon_cks[1, :, :, :] - ucon_cks[0, :, :, :] * gcon_cks[0, 1, :, :, :] / gcon_cks[0, 0, :, :, :]
    Uprims[:, :, :, 1] = ucon_cks[2, :, :, :] - ucon_cks[0, :, :, :] * gcon_cks[0, 2, :, :, :] / gcon_cks[0, 0, :, :, :]
    Uprims[:, :, :, 2] = ucon_cks[3, :, :, :] - ucon_cks[0, :, :, :] * gcon_cks[0, 3, :, :, :] / gcon_cks[0, 0, :, :, :]

    print(" - reconstructing B primitives")
    Bprims = np.zeros((n1, n2, n3, 3))
    Bprims[:, :, :, 0] = np.einsum('abc,abc->abc',bcon_cks[1],ucon_cks[0])-np.einsum('abc,abc->abc',bcon_cks[0],ucon_cks[1])
    Bprims[:, :, :, 1] = np.einsum('abc,abc->abc',bcon_cks[2],ucon_cks[0])-np.einsum('abc,abc->abc',bcon_cks[0],ucon_cks[2])
    Bprims[:, :, :, 2] = np.einsum('abc,abc->abc',bcon_cks[3],ucon_cks[0])-np.einsum('abc,abc->abc',bcon_cks[0],ucon_cks[3])

    ## begin remapping
    print(f" - remapping to {target_n1}x{target_n2}x{target_n3}")

    extrema = np.abs(np.array([all_x1s.min(), all_x1s.max(), all_x2s.min(), all_x2s.max(), all_x3s.min(), all_x3s.max()]))
    r_min = reh * 0.95
    r_max = extrema.min() * 0.98
    r_lin = np.logspace(np.log10(r_min), np.log10(r_max), target_n1)
    h_lin = np.linspace(0, np.pi, target_n2+1)
    h_lin = (h_lin[1:] + h_lin[:-1]) / 2.
    p_lin = np.linspace(0, 2.*np.pi, target_n3+1)[:-1]
    R_ks, H_ks, P_ks = np.meshgrid(r_lin, h_lin, p_lin, indexing='ij')

    # get CKS locations for R_ks, H_ks, P_ks
    X_cks = R_ks * np.cos(P_ks) * np.sin(H_ks) - bhspin * np.sin(P_ks) * np.sin(H_ks)
    Y_cks = R_ks * np.sin(P_ks) * np.sin(H_ks) + bhspin * np.cos(P_ks) * np.sin(H_ks)
    Z_cks = R_ks * np.cos(H_ks)

    ## interpolate the CKS U & B
    print(" - interpolating primitive U")
    Uprims_sph = np.zeros((target_n1, target_n2, target_n3, 3))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Uprims[:,:,:,0].transpose((2, 1, 0)), method=interp_method)
    Uprims_sph[:, :, :, 0] = rgi((X_cks, Y_cks, Z_cks))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Uprims[:,:,:,1].transpose((2, 1, 0)), method=interp_method)
    Uprims_sph[:, :, :, 1] = rgi((X_cks, Y_cks, Z_cks))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Uprims[:,:,:,2].transpose((2, 1, 0)), method=interp_method)
    Uprims_sph[:, :, :, 2] = rgi((X_cks, Y_cks, Z_cks))

    print(" - interpolating primitive B")
    Bprims_sph = np.zeros((target_n1, target_n2, target_n3, 3))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Bprims[:,:,:,0].transpose((2, 1, 0)), method=interp_method)
    Bprims_sph[:, :, :, 0] = rgi((X_cks, Y_cks, Z_cks))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Bprims[:,:,:,1].transpose((2, 1, 0)), method=interp_method)
    Bprims_sph[:, :, :, 1] = rgi((X_cks, Y_cks, Z_cks))
    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), Bprims[:,:,:,2].transpose((2, 1, 0)), method=interp_method)
    Bprims_sph[:, :, :, 2] = rgi((X_cks, Y_cks, Z_cks))

    ## reconstruct four-vectors in CKS
    print(" - reconstructing ucon in CKS at new grid")
    gcon_cks_interp = metrics.cks_inverse_metric(X_cks, Y_cks, Z_cks, bhspin)
    gcov_cks_interp = metrics.cks_metric(X_cks, Y_cks, Z_cks, bhspin)
    alpha = 1. / np.sqrt(-gcon_cks_interp[0,0,:,:,:])

    gamma = np.sqrt(1. + np.einsum('abci,abci->abc',
                                   np.einsum('ijabc,abci->abcj',gcov_cks_interp[1:,1:,:,:,:],Uprims_sph),
                                   Uprims_sph))
    ucon_cks = np.zeros((target_n1, target_n2, target_n3, 4))
    ucon_cks[:,:,:,1:] = Uprims_sph - gamma[:,:,:,None] * alpha[:,:,:,None] * gcov_cks_interp[0,1:,:,:,:].transpose((1,2,3,0))
    ucon_cks[:,:,:,0] = gamma / alpha
    ucov_cks = np.einsum('ijabc,abci->abcj', gcov_cks_interp, ucon_cks)
    # usq = np.einsum('abci,abci->abc', ucon_cks, ucov_cks)  # should be ~ -1

    print(" - reconstructing bcon")
    bcon_cks = np.zeros((target_n1, target_n2, target_n3, 4))
    bcon_cks[:,:,:,0] = np.einsum('abci,abci->abc', Bprims_sph, ucov_cks[:,:,:,1:])
    bcon_cks[:,:,:,1:] = ( Bprims_sph + ucon_cks[:,:,:,1:] * bcon_cks[:,:,:,0,None] ) / ucon_cks[:,:,:,0,None]
    bcov_cks = np.einsum('ijabc,abci->abcj', gcov_cks_interp, bcon_cks)
    # bdotu = np.einsum('abci,abci->abc', bcon_cks, ucov_cks)  # should be ~ 0

    ## then translate to eks via ks
    print(" - translating ucon and bcon to EKS via KS")
    ucon_ks = metrics.cks_vec_to_ks(ucon_cks.transpose((3,0,1,2)), X_cks, Y_cks, Z_cks, a=bhspin).transpose((1,2,3,0))
    bcon_ks = metrics.cks_vec_to_ks(bcon_cks.transpose((3,0,1,2)), X_cks, Y_cks, Z_cks, a=bhspin).transpose((1,2,3,0))
    dXdx = np.zeros((target_n1, target_n2, target_n3, 4, 4))
    dXdx[:, :, :, 0, 0] = 1.
    dXdx[:, :, :, 1, 1] = 1. / R_ks
    dXdx[:, :, :, 2, 2] = 1.
    dXdx[:, :, :, 3, 3] = 1.
    ucon_eks = np.einsum('abciu,abci->abcu', dXdx, ucon_ks)
    bcon_eks = np.einsum('abciu,abci->abcu', dXdx, bcon_ks)

    ## now pack up primitives
    prims = np.zeros((target_n1, target_n2, target_n3, 8))

    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), rho.transpose((2, 1, 0)), method=interp_method)
    prims[:, :, :, 0] = rgi((X_cks, Y_cks, Z_cks))

    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), internal_energy.transpose((2, 1, 0)), method=interp_method)
    prims[:, :, :, 1] = rgi((X_cks, Y_cks, Z_cks))

    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), pressure_mag.transpose((2, 1, 0)), method=interp_method)
    new_pressure_mag = rgi((X_cks, Y_cks, Z_cks))

    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), beta.transpose((2, 1, 0)), method=interp_method)
    new_beta = rgi((X_cks, Y_cks, Z_cks))

    rgi = RegularGridInterpolator((all_x1s, all_x2s, all_x3s), sigma.transpose((2, 1, 0)), method=interp_method)
    new_sigma = rgi((X_cks, Y_cks, Z_cks))

    gcon_eks = metrics.get_gcon_eks_3d(R_ks, H_ks, a=bhspin)
    prims[:, :, :, 2] = ucon_eks[:, :, :, 1] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 1] / gcon_eks[:, :, :, 0, 0]
    prims[:, :, :, 3] = ucon_eks[:, :, :, 2] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 2] / gcon_eks[:, :, :, 0, 0]
    prims[:, :, :, 4] = ucon_eks[:, :, :, 3] - ucon_eks[:, :, :, 0] * gcon_eks[:, :, :, 0, 3] / gcon_eks[:, :, :, 0, 0]

    prims[:, :, :, 5] = bcon_eks[:, :, :, 1] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 1]
    prims[:, :, :, 6] = bcon_eks[:, :, :, 2] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 2]
    prims[:, :, :, 7] = bcon_eks[:, :, :, 3] * ucon_eks[:, :, :, 0] - bcon_eks[:, :, :, 0] * ucon_eks[:, :, :, 3]

    ## output
    ofname = fname.replace(".athdf", "") + ".h5"
    ohfp = h5py.File(ofname, 'w')
    write_header(ohfp, R_ks, H_ks, P_ks, fluid_gamma, bhspin)
    ohfp['prims'] = prims
    ohfp['beta'] = new_beta
    ohfp['sigma'] = new_sigma
    ohfp['pressure_mag'] = new_pressure_mag
    ohfp.close()


