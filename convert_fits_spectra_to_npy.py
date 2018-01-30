from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import glob

home = os.getenv('HOME')
stacking_analysis_dir = home + "/Desktop/FIGS/stacking-analysis-pears/"

sys.path.append(stacking_analysis_dir + 'codes/')
import fast_chi2_jackknife as fcj

if __name__ == '__main__':
    
    # define directories
    models_dir = home + '/Documents/GALAXEV_BC03/bc03/src/cspout_new/'
    metallicity = 'm62/'

    for fl in glob.glob(models_dir + metallicity + '*.fits'):

        # get just the file name
        fl_name = os.path.basename(fl)
        print "Working on", fl_name

        # open file and save all data to npy arrays
        # for each fits file there will be three numpy arrays
        # one npy array for the lambda grid
        # one npy array for all ages
        # one npy array for all spectra which are sampled on the corresponding 
        # lambda grid and there is one spectrum for each age
        hdulist = fits.open(fl)
        
        # save lambda grid and ages separately
        npy_lamgrid_to_save = hdulist[1].data
        npy_ages_to_save = hdulist[2].data

        np.save(fl.replace('.fits', '_lamgrid.npy'), npy_lamgrid_to_save)
        np.save(fl.replace('.fits', '_ages.npy'), npy_ages_to_save)

        # now loop over all spectra and save them in a single numpy array
        total_ext = fcj.get_total_extensions(hdulist)
        npy_spectra_to_save = np.empty([total_ext-2, len(npy_lamgrid_to_save)])
        for i in range(total_ext - 2):
            npy_spectra_to_save[i] = hdulist[i+3].data

        np.save(fl.replace('.fits', '_allspectra.npy'), npy_spectra_to_save)

        # close fits file
        hdulist.close()

    sys.exit(0)