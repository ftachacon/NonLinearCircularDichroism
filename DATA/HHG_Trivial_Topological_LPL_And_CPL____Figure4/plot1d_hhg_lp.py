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

htrivial = np.loadtxt("HHG__Phi0__0.060__M0t2__2.540__e0__0.000__E00.0045__.dat")

htopo = np.loadtxt("HHG__Phi0__1.160__M0t2__2.540__e0__0.000__E00.0045__.dat")

Ntime           = int(htrivial[0,0]);
dt              = htrivial[0,5];

print( "Ntime = ", Ntime, ";  dt = ",dt);


om      = htrivial[1:Ntime,0];
hhgtriv = htrivial[1:Ntime,5];
hhgtopo = htopo[1:Ntime,5];


#################################
#################################
width   = 11.4;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.20, 0.162, 0.77, 0.80]);

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);

hzero   = 0.25e-16;

p1,=plt.plot( om, hhgtriv + hzero, 'b', lw=2   );
p2,=plt.plot( om, hhgtopo + hzero, '-r', lw=1.1 );

plt.fill_between( om, hhgtopo, hhgtriv, color=[0.4,0.,0.95], alpha=0.250 )


#plt.legend([p1, p2], [ r'$C \,=\,\,\,\, 0$', r'$C \,=\, +1$'],fontsize=28);

plt.legend([p1, p2], [r'$\rm \,\,Trivial\,\,\,\,$', r'$\rm Topological$'],fontsize=31);
xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

plt.text(.32, 5e-2, '(a) Linear Driver', fontsize=34, color="k");
#plt.text(2, 5e-2, '      Linear Driver', fontsize=28, color="k");




plt.yscale('log')
plt.xlim([0,33])

plt.grid(True)

plt.xlabel(r'$\rm Harmonic\,\,Order$', fontsize=37 );
#plt.ylabel(r'$\rm I_{HHG} (arb.\, u.)$', fontsize=34 );

plt.text(-.256, 0.5, r'$\rm I_{HHG}\,\,(arb.\,u.)$ ',
         fontsize=37,
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes)


plt.tick_params( labelsize = 37 );

yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.0e0])


#fname='HHG_Spectrum_Trivial_Vs_NonTrivial__Linear__Polarization.pdf'

fname='Figure4a.pdf'
plt.savefig(fname, dpi = 300);


plt.show()



