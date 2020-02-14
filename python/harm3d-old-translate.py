# Translate 

import os
import h5py
import numpy as np
import sys

if __name__ == "__main__":

    fnames = sys.argv[1:]
    
    for fname in fnames:
        print("translating {0:s}".format(fname))
        # old harm3d output is HDF5!
        ofname = os.path.join("dumps", fname)
        
        # parse header
        fp = h5py.File(fname, 'r')
        
        grid = fp['Header/Grid']
        # Grid variables
        # Old N1,2,3 are for a single process
        N1 = grid['totalsize1'][0]
        N2 = grid['totalsize2'][0]
        N3 = grid['totalsize3'][0]
        NG = grid['NG'][0]
        dx1 = grid['dx1'][0]
        dx2 = grid['dx2'][0]
        dx3 = grid['dx3'][0]
        # Hotaka uses startx for the *ghost* zones.
        # But doesn't output them
        # So we just change the convention here
        startx1 = grid['startx1'][0] + NG * dx1
        startx2 = grid['startx2'][0] + NG * dx2
        startx3 = grid['startx3'][0] + NG * dx3
        
        # Domain variables
        t = grid['t'][0]
        a = grid['a'][0]
        gam = grid['gam'][0]
        rin = grid['Rin'][0]
        rout = grid['Rout'][0]
        # Modern IL convention is hslope=1 corresponds to th=x2*pi
        hslope = 1 - grid['h_slope'][0]
        R0 = grid['R0'][0]
        
        # IO variables
        hdr_io = fp['Header/IO']
        dump_cadence = hdr_io['DT_dump_out'][0]
        n_dump = hdr_io['N_out'][0]
        n_step = hdr_io['nstep'][0]
        
        # Get prims data, while verifying some of the above is true
        prims = np.zeros((N1, N2, N3, 8))
        prims[:,:,:,0] = fp['rho'][()]
        prims[:,:,:,1] = fp['uu'][()]
        prims[:,:,:,2] = fp['v1'][()]
        prims[:,:,:,3] = fp['v2'][()]
        prims[:,:,:,4] = fp['v3'][()]
        prims[:,:,:,5] = fp['B1'][()]
        prims[:,:,:,6] = fp['B2'][()]
        prims[:,:,:,7] = fp['B3'][()]
        
        # get output file
        hfp = h5py.File(ofname, 'w')
        
        # write header
        hfp.create_group("/header")
        hfp.create_group("/header/geom")
        hfp.create_group("/header/geom/mks")
        hfp.create_group("/header/units")
        
        hdr = hfp['header']
        hdr['has_electrons'] = 0
        
        hdr['n1'] = N1
        hdr['n2'] = N2
        hdr['n3'] = N3
        hdr['n_prim'] = 8
        hdr['gam'] = gam
        
        geom = hdr['geom']
        geom['startx1'] = startx1
        geom['startx2'] = startx2
        geom['startx3'] = startx3
        geom['dx1'] = dx1
        geom['dx2'] = dx2
        geom['dx3'] = dx3
        
        # TODO does this == b"MKS"?  Write strings from python good
        hfp.create_dataset("/header/metric", data=np.string_("MKS"))
        geom['mks']['a'] = a
        geom['mks']['r_in'] = rin
        geom['mks']['r_out'] = rout
        geom['mks']['R0'] = R0
        geom['mks']['hslope'] = hslope
        
        # Per-dump numbers
        hfp['dump_cadence'] = dump_cadence
        hfp['full_dump_cadence'] = -1
        hfp['t'] = t
        
        hfp['is_full_dump'] = 0
        hfp['n_dump'] = n_dump
        hfp['n_step'] = n_step
        
        hfp['prims'] = prims
        hfp.close()
