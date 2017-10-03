from __future__ import division
import numpy as np
import glob, os

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
pgf_preamble = {"pgf.texsystem": "pdflatex"}
mpl.rcParams.update(pgf_preamble)


if __name__ == '__main__':
    cspout = '/Users/baj/Documents/GALAXEV_BC03/bc03/src/cspout_new/'
    metals = ["m22","m32","m42","m52","m62","m72"]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plotfiles = []
    reds = plt.cm.get_cmap('Reds')
    
    for metallicity in metals:
        metalsfolder = cspout + metallicity + '/'
        for file in glob.glob(metalsfolder + "*_tauV10_*_tau5_*.3color"):
            #plotfiles.append(file)
            
            if file.split('/')[-1].split('_')[2] == 'm22': Z = 0.0001
            if file.split('/')[-1].split('_')[2] == 'm32': Z = 0.0004
            if file.split('/')[-1].split('_')[2] == 'm42': Z = 0.004
            if file.split('/')[-1].split('_')[2] == 'm52': Z = 0.008
            if file.split('/')[-1].split('_')[2] == 'm62': Z = 0.02
            if file.split('/')[-1].split('_')[2] == 'm72': Z = 0.05

            breakampfile = np.genfromtxt(file, dtype =None, skip_header=29, names=['age', 'dn4000'], usecols=(0,2))
            physpropfile = np.genfromtxt(file.replace('3color', '4color'), dtype =None, skip_header=29,\
                                         names=['age', 'bmag', 'vmag', 'sfr'], usecols=(0,2,3,9))

            filename = os.path.basename(file)
            tauV = float(filename.split('_')[3][4:]) / 10
            dn4000 = breakampfile['dn4000']

            print filename, tauV

            plt.scatter(Z*np.ones(len(dn4000)), breakampfile['age'], marker='o', s=75,\
                        c=dn4000, vmin=min(dn4000), vmax=max(dn4000), cmap=reds, edgecolor='')

#plt.contour(Z*np.ones(len(dn4000)), breakampfile['age'], dn4000)

    plt.colorbar()
    plt.show()