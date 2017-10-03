from __future__ import division
import numpy as np
import glob

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

if __name__ == '__main__':
    cspout = '/Users/baj/Documents/GALAXEV_BC03/bc03/src/cspout_new/'
    #metals = ["m22","m32","m42","m52","m62","m72"]
    metals = ["m32"]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plotfiles = []
    ur_col = []
    stellarmass = []
    
    for metallicity in metals:
        metalsfolder = cspout + metallicity + '/'
        for file in glob.glob(metalsfolder + "*.1color"):
            if "tau0" in file: continue
            
            ubv_col = np.genfromtxt(file, dtype =None, skip_header=29, names=['age', 'umag', 'bmag', 'vmag'], usecols=(0,2,3,4))
            rmag = np.genfromtxt(file.replace('1color','2color'), dtype =None, skip_header=29, names=['rm'], usecols=(1))
            mstar = np.genfromtxt(file.replace('1color','4color'), dtype =None, skip_header=29, names=['ms'], usecols=(6))

            for i in range(len(ubv_col)):
                ur = ubv_col['umag'][i] - rmag['rm'][i]
                ur_col.append(ur)

                stellarmass.append(mstar['ms'][i])

    stellarmass = np.array(stellarmass)
    ur_color = np.array(ur_col)

    ax.set_xlabel(r'$\frac{M_*}{[1\ M_\odot]}$')
    ax.set_ylabel('$\mathrm{(U-R)_{rest}}$')
    ax.plot(stellarmass, ur_color, 'o', markersize=2, color='k', markeredgecolor='none')

    ax.set_ylim(-2.0, 3.0)
    ax.set_xlim(0.0, 0.8)

    plt.show()