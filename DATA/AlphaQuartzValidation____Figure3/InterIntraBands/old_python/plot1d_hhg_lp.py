#!/usr/bin/python
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import numpy as np
import os
import sys
import shutil



import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker

from matplotlib.pylab import *
from numpy import arange
from scipy.interpolate import interp2d


plt.rc('lines', lw=3, color='k')

htrivial = np.loadtxt("HHG__Phi0__0.000__M0t2__12.375__e0__0.000__.dat")


Ntime           = int(htrivial[0,0]);
dt              =   htrivial[0,5];

print "Ntime = ", Ntime, ";  dt = ",dt;


om      = htrivial[1:Ntime,0];
hhgtriv = htrivial[1:Ntime,5];


#################################
#################################
width   = 11;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2.0);

hzero   = 0.995e-16;

p1,=plt.plot( om, hhgtriv, 'r', lw=2   );
plt.fill_between( om, hhgtriv, 1.e-16, color=[0.8,0,0], alpha=0.25 )


#p2,=plt.plot( om, hhgtopo + hzero, 'r', lw=1.1 );
plt.legend([p1], [r'$\alpha$-'+'quartz'] ,fontsize=26);
xticks0  = np.arange(3,50,4);
plt.xticks( xticks0 );

plt.text(3.5, 1.5e-2, '(b)', fontsize=34, color="k");
plt.text(3.5, 1.5e-2, '        Theory', fontsize=28, color="k");

plt.yscale('log')
plt.xlim([3.,25])

plt.grid(True)

#plt.xlabel(r'$\rm Harmonic-Order$', fontsize=34 );
#plt.ylabel(r'$\rm I_{HHG} (arb.\, u.)$', fontsize=34 );
plt.tick_params(labelsize = 30 );

yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.00e-0])

#fname='HHG_Spectrum_Trivial_Vs_NonTrivial.pdf';
#filename1           = FigureDir + fname;
#fileNamePicture     = ProjPath  + filename1;
#plt.savefig(fname, dpi = 300);


fname='Theory__For__Woerner__Experiment.pdf'
#filename1           = FigureDir + fname;
#fileNamePicture     = ProjPath  + filename1;
plt.savefig(fname, dpi = 400);


plt.show()



