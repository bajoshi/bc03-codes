"""
    This code generates plots that show the evolution of D_n(4000) with age for a CSP.
    Metallicity is fixed for each run and the spectra are shown within a range of tau values
    i.e. 0.01 <= tau (Gyr) <= 10.01 in steps of 50 x 10^6 years (0.05 Gyr).
"""

from __future__ import division
import numpy as np
import array
import pyfits as pf
import glob

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
pgf_preamble = {"pgf.texsystem": "pdflatex"}
mpl.rcParams.update(pgf_preamble)

# Cosmology used
H_0 = 69.6
omega_m0 = 0.286
omega_r0 = 8.24e-5
omega_lam0 = 0.714

# other CONSTANTS
solar_lum = 3.826e33 # erg s^-1
c = 3e18 # A s^-1

if __name__ == '__main__':
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('$\mathrm{log}\,(t/\mathrm{yr})$', fontsize=18)
    ax1.set_ylabel('$\mathrm{log}\,(\mathrm{D_n}(4000))$', fontsize=18)

    cspout = '/Users/baj/Documents/GALAXEV_BC03/bc03/src/cspout/'
    metalsfolder = 'm62/'
    files = glob.glob(cspout + metalsfolder + '*1_salp.fits')
    files = [file.split('/')[-1] for file in files]
    sortedfiles = sorted(files, key=lambda name: float(name[21:-10]))

    tauarr = np.arange(0.01, 10.05, 0.1)

    cubehelix = cm.cubehelix(np.arange(len(tauarr)))

    count = 0
    for file in sortedfiles:

        tau = float(file.split('/')[-1].split('_')[4][3:])
        print tau

        hdu = pf.open(cspout + metalsfolder + file, memmap=False)
    
        rest_lam = hdu[1].data
        agevals = hdu[2].data
        
        arg4100 = np.argmin(abs(rest_lam - 4100))
        arg4000 = np.argmin(abs(rest_lam - 4000))
        arg3950 = np.argmin(abs(rest_lam - 3950))
        arg3850 = np.argmin(abs(rest_lam - 3850))

        break_measured = []

        for i in range(3,248,1):
            spec = hdu[i].data
        
            lum = spec * solar_lum
            
            sum_up = np.trapz(lum[arg4000:arg4100], x=rest_lam[arg4000:arg4100])
            sum_low = np.trapz(lum[arg3850:arg3950], x=rest_lam[arg3850:arg3950])
            
            break_measured.append(sum_up/sum_low)

        ax1.plot(np.log10(agevals), np.log10(break_measured), color=cubehelix[count], label=str(tau))
        count += 1
        hdu.close(verbose=True,closed=True)

    ax1.legend(loc=0, numpoints=1, prop={'size':8}, ncol=5)
    # size 5 and ncol 7 seems to work for most

    """ 
        Fix the limits on Y axis -- it should be the same for all plots to compare.
        Change it to -0.05 to 0.40
    """

    ax1.axvline(x=np.log10(7.948e9), color='k', linestyle='--')
    ax1.axvline(x=np.log10(13.8e9), color='k', linestyle='--')
    ax1.axhline(y=0, linestyle='--')
    fig1.savefig('break_m62_tau', dpi=300)
    #plt.show()
