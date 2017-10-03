"""
    This code is largely based on the code for reading ised files 
    from the ezgal package. 
    See the ezgal github page
    and also refer to Mancone & Gonzalez 2012, PASP, 124, 606.
"""

from __future__ import division

import numpy as np
import array
from astropy.io import fits as pf

import glob 
import os 
import sys


def read_current_filepos(filehandle, type='i', number=1):
    """
        This function reads the given number of binary data which is of the supplied type
        from the specified filehandle. The default number of elements to extract is one
        and the default type is an integer.
    """

    arr = array.array(type)
    arr.fromfile(filehandle, number)

    return np.asarray(arr)

def read_ised(modelfile, del_modelfile):
    
    fh = open(modelfile, 'rb')
    
    # For this file
    # int = 4 bytes
    # float = 4 bytes
    # char = 1 byte

    ignore = read_current_filepos(fh) # this is the number 1208 for SSPs and depends on Tcut for CSPs
    totalages = read_current_filepos(fh)[0] # this should be 221 for the SSP models
    
    """
        FOR CSPs ----
        The number of ages you get in the previous line depends on the CUT OFF TIME set in the csp run.
        The max age of any model that is given by BC03 is 20 Gyr. I found by trial and error that if the
        ---> Tcut is set to above 20 Gyr then it gives 221 ages which is the same as that given by a SSP.
        ---> Tcut is set to undet 20 Gyr but more than 13.8 Gyr then it gives a number (for totalages)
             that is between 221 and 245.
        ---> Tcut is set to 13.8 Gyr then it gives 245 ages.
        ---> Tcut is set to 5 Gyr then it gives 261 ages.
        Keep in mind however that the total age range still remains between 0 to 20 Gyr.
        Also I'm not sure if having a lower/higher number of totalages (than 245) will have any effect on my analysis.
        For now I have kept Tcut at 13.8 Gyr.
    """
    
    allages = read_current_filepos(fh, type='f', number=totalages) # in years
    
    fh.seek(328, 1) # done by trial and error and looking at the ezgal code but mostly the ezgal code helped
    # this goes 328 bytes forward from the current position
    # if the second keyword is not provided then it assumes absolute positioning

    totalwavelengths = read_current_filepos(fh)[0]
    allwavelengths = read_current_filepos(fh, type='f', number=totalwavelengths)
    
    seds = np.zeros((totalages, totalwavelengths))
    
    hdu = pf.PrimaryHDU()
    
    tau = int(modelfile.split('/')[-1].split('_')[5][3:])
    tauV = int(modelfile.split('/')[-1].split('_')[3][4:])
    
    hdulist = pf.HDUList(hdu)
    hdulist.append(pf.ImageHDU(allwavelengths))
    hdulist.append(pf.ImageHDU(allages))
    
    for i in range(totalages):
        ignore = read_current_filepos(fh, number=2)
        nlam = read_current_filepos(fh)
    
        seds[i] = read_current_filepos(fh, type='f', number=totalwavelengths)
        
        no = read_current_filepos(fh)
        ignore = read_current_filepos(fh, type='f', number=no)
    
        hdulist.append(pf.ImageHDU(seds[i]))
    
    hdulist.writeto(modelfile.replace('.ised','.fits'), clobber=True)
    print "Writing ...", modelfile.replace('.ised','.fits')
    
    if del_modelfile:
        os.remove(modelfile)
        print "Deleted ...", modelfile

    fh.close()

if __name__ == '__main__':
    """
    # for SSPs
    for file in glob.glob('/Users/baj/Documents/GALAXEV_BC03/bc03/models/Padova1994/salpeter/' + '*.ised'):
        read_ised(file)
    """
    # for CSPs
    cspout = '/Users/baj/Documents/GALAXEV_BC03/bc03/src/cspout_new/'
    metals = ["m22", "m32","m42","m52","m62","m72"]
    
    for metallicity in metals:
        metalsfolder = metallicity + '/'
        dir = cspout + metalsfolder
        for file in glob.glob(dir + '*csp*.ised'):
            read_ised(file)
