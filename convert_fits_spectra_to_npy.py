from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import glob

home = os.getenv('HOME')

if __name__ == '__main__':
    
    # define directories
    models_dir = home + '/Documents/GALAXEV_BC03/bc03/src/cspout_new/'
    metallicity = 'm62/'

    for fl in glob.glob(models_dir + metallicity + '*.fits'):

        fl_name = os.path.basename(fl)
        

    sys.exit(0)