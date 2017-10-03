from __future__ import division

import numpy as np
import array

from astropy.io import fits as pf

import glob
import sys
import os

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end

import astropy.units as u
from astropy.cosmology import z_at_value
from astropy.cosmology import Planck15

if __name__ == '__main__':

    # Cosmology used
    H_0 = 69.6
    omega_m0 = 0.286
    omega_r0 = 8.24e-5
    omega_lam0 = 0.714
    
    # other CONSTANTS
    solar_lum = 3.826e33 # erg s^-1
    c = 3e18 # A s^-1
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel(r'$\mathrm{log}\,(t/\mathrm{yr})$', fontsize=18)
    ax1.set_ylabel(r'$\mathrm{D_n}(4000)$', fontsize=18)
    
    if __name__ == '__main__':
        for file in glob.glob(home + '/Documents/GALAXEV_BC03/bc03/models/Padova1994/salpeter/' + '*.fits'):
            
            if 'm22' in file:
                metals = 0.0001
                metallabel = r'$0.005\, \mathrm{Z}_\odot$'
            elif 'm32' in file:
                metals = 0.0004
                metallabel = r'$0.02\, \mathrm{Z}_\odot$'
            elif 'm42' in file:
                metals = 0.004
                metallabel = r'$0.2\, \mathrm{Z}_\odot$'
            elif 'm52' in file:
                metals = 0.008
                metallabel = r'$0.4\, \mathrm{Z}_\odot$'
            elif 'm62' in file:
                metals = 0.02
                metallabel = r'$\mathrm{Z}_\odot$'
            elif 'm72' in file:
                metals = 0.05
                metallabel = r'$2.5\, \mathrm{Z}_\odot$'
            
            hdu = pf.open(file)
    
            rest_lam = hdu[1].data
            agevals = hdu[2].data
    
            arg4100 = np.argmin(abs(rest_lam - 4100))
            arg4000 = np.argmin(abs(rest_lam - 4000))
            arg3950 = np.argmin(abs(rest_lam - 3950))
            arg3850 = np.argmin(abs(rest_lam - 3850))
    
            break_measured = []
    
            for i in range(3,224,1):
                spec = hdu[i].data
                
                lum = spec * solar_lum
    
                sum_up = np.trapz(lum[arg4000:arg4100], x=rest_lam[arg4000:arg4100])
                sum_low = np.trapz(lum[arg3850:arg3950], x=rest_lam[arg3850:arg3950])
    
                break_measured.append(sum_up/sum_low)

            break_measured = np.asarray(break_measured)
            ax1.plot(np.log10(agevals), break_measured, label=metallabel)

            age_indices = np.where((agevals>=3e9) & (agevals<=9e9))[0]

            if metals == 0.02:
                solar_val_break = break_measured[age_indices]
            if metals == 0.05:
                highest_metal_val_break = break_measured[age_indices]

    # run this part if you want the maximum break value (assuming solar mettalicity SSP) at a given redshift
    # get maximum dn4000 at age
    ages = np.linspace(3,9,len(solar_val_break))*u.Gyr
    redshiftvals = [z_at_value(Planck15.age, age) for age in ages]

    for i in range(len(highest_metal_val_break)):
        #print solar_val_break[i], redshiftvals[i]
        print highest_metal_val_break[i], redshiftvals[i]

    #fig = plt.figure()
    #ax = fig.add_subplot(111)

    #ax.plot(np.log10(ages.value*1e9), solar_val_break)
    #plt.show()
 
    ax1.minorticks_on()
    ax1.tick_params('both', width=1, length=3, which='minor')
    ax1.tick_params('both', width=1, length=4.7, which='major')
    
    ax1.set_xlim(6.5,np.log10(13.72e9))
    
    ax1.legend(loc=0, frameon=False)
    #fig1.savefig('break_ssp.eps', dpi=300)
    
    #plt.show()
    sys.exit(0)


