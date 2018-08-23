"""
 This code is based on the code run_grid.py
 very similar but this one is to be used to 
 extract only a sinlge spectrum.
 Hopefully, this is a little better commented as well.
 
 This code HAS to be kept in the same folder as the ised files. 
 This is because the actual bc03 code written in fortran is in there
 and csp_galaxev needs to be able to use these codes.
 
 None of the bc03 programs (and this code which calls the bc03
 programs) will work if you have not used >> make
 to compile the programs on the computer being used.

 0. This code also HAS to be run in csh.
 
 1. If you run the code and nothing happens then you need to 
 source .bc_cshrc again. You will have to source .bc_cshrc 
 everytime this happens.

 2. If it gives the error 

 Please assign correct filter file
 Use command: stdfilt

 before even beginning to do anything in the program then you will 
 have to >>source .bc_cshrc again.

"""

from __future__ import division

import numpy as np
from astropy.io import fits

import subprocess
import os
import sys
import shutil

home = os.getenv('HOME')
bc03_basedir = home + '/Documents/GALAXEV_BC03/'

sys.path.append(bc03_basedir + 'bc03-codes/')
import bc03_extractspectra as bce

def callcsp(isedfile, metallicity, tauV, mu, tau, dust='n', sfh='1', recycle='y', ep='0.1', tcut='13.8'):
    """
    This function expects to get tau in units of Gyr.
    Set the rest of the parameters below.
    """

    start = subprocess.Popen([home + '/Documents/GALAXEV_BC03/bc03/src/csp_galaxev'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    output = "bc2003_hr_" + metallicity + "_tauV" + str(int(tauV*10)) + "_csp_tau" + str(int(float(str(tau)[0:6])*10000)) + "_salp"
    print output
    params = start.communicate(os.linesep.join([isedfile, dust, str(tauV), str(mu), sfh, str(tau), recycle, ep, tcut, output]))
    start.poll()

def removefiles(files):
    for f in files:
        os.remove(f)

if __name__ == '__main__':

    # Set your parameters here
    metals = 0.0001  # absolute fraction of metals
    tau = 5  # in gyr
    tauV = 2.0
    mu = 0.3
    savedir = bc03_basedir
    # if you want to extract a single age then enter the age below
    # it will extract the spectrum for this age and save it to a text file.
    age = 5e9  # in years

    if metals == 0.0001:
        metallicity = 'm22'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m22_salp_ssp.ised"
    elif metals == 0.0004:
        metallicity = 'm32'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m32_salp_ssp.ised"
    elif metals == 0.004:
        metallicity = 'm42'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m42_salp_ssp.ised"
    elif metals == 0.008:
        metallicity = 'm52'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m52_salp_ssp.ised"
    elif metals == 0.02:
        metallicity = 'm62'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m62_salp_ssp.ised"
    elif metals == 0.05:
        metallicity = 'm72'
        isedfile = home + "/Documents/GALAXEV_BC03/bc03/src/bc2003_hr_m72_salp_ssp.ised"

    callcsp(isedfile, metallicity, tauV, mu, tau, dust='n', sfh='1', recycle='y', ep='0.1', tcut='13.8')

    # convert ised file to fits
    isedname = "bc2003_hr_" + metallicity + "_tauV" + str(int(tauV*10)) + "_csp_tau" + str(int(float(str(tau)[0:6])*10000)) + "_salp.ised"
    bce.read_ised(isedname, del_modelfile=False)

    # remove files that aren't required 
    isedbasename = isedname.split('.')[0]
    files = []
    files.append(isedbasename + '.1ABmag')
    files.append(isedbasename + '.1color')
    files.append(isedbasename + '.2color')
    files.append(isedbasename + '.3color')
    files.append(isedbasename + '.4color')
    files.append(isedbasename + '.5color')
    files.append(isedbasename + '.6lsindx_ffn')
    files.append(isedbasename + '.6lsindx_sed')
    files.append(isedbasename + '.6lsindx_sed_lick_system')
    files.append(isedbasename + '.7lsindx_ffn')
    files.append(isedbasename + '.7lsindx_sed')
    files.append(isedbasename + '.7lsindx_sed_lick_system')
    files.append(isedbasename + '.8lsindx_sed_fluxes')

    removefiles(files)

    # move the ised and fits files
    shutil.move(bc03_basedir + 'bc03/src/' + isedname, savedir + isedname)
    ised_fitsname = isedname.replace('.ised', '.fits')
    shutil.move(bc03_basedir + 'bc03/src/' + ised_fitsname, savedir + ised_fitsname)

    # extract single spectrum
    hdulist = fits.open(savedir + ised_fitsname)
    ages = hdulist[2].data
    wav = hdulist[1].data

    age_idx = np.argmin(abs(ages - age))
    spec = hdulist[age_idx + 3].data

    # save it in a txt file
    data = np.array(zip(spec, wav), dtype=[('spec', float), ('wav', float)])
    np.savetxt(savedir + 'spec_age' + str('{:.2e}'.format(age)).replace('.', 'p') + '_metals' + str(metals).replace('.', 'p') + '.txt', \
        data, fmt=['%.4e', '%.4f'], delimiter=' ', header='tau=' + str(tau) + ' Gyr, tau_V=' + str(tauV))

    sys.exit(0)