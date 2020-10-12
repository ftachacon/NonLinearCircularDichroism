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

htrivial = np.loadtxt("HHG__Phi0__0.060__M0t2__2.540__e0__-1.000__.dat");
htopo = np.loadtxt("HHG__Phi0__1.160__M0t2__2.540__e0__-1.000__.dat");

Ntime           = int(htrivial[0,0]);
dt              =   htrivial[0,5];

print ("Ntime = ", Ntime, ";  dt = ",dt);


om          = htrivial[1:Ntime,0];
hhgtriv     = htrivial[1:Ntime,5];
hhgtopo     = htopo[1:Ntime,5];




#################################
#################################
width   = 11.3;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.15, 0.162, 0.81, 0.81]);

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);

hzero   = 0.25e-16;

p1,=plt.plot( om, hhgtriv + hzero, 'b', lw=2   );
p2,=plt.plot( om, hhgtopo + hzero, '-r', lw=2 );

plt.fill_between( om, hhgtopo, hhgtriv, color=[1,1,1], alpha=0.71 )

h1 = 0.50e-6
h2 = 0.5e-10
#plt.plot([0,16.25],[h1,h1],color=[0.85,0.,0],lw=4,linestyle='--')
#plt.plot([0,19.2],[h2,h2],color=[0.,0.,0.95],lw=4,linestyle='--')

#plt.annotate('Top. cutoff', xy=(16, 0.5e-6), fontsize=31,color=[0.8,0.,0.],  xycoords='data',
#            xytext=(0.6, 0.7), textcoords='axes fraction',
#            arrowprops=dict(facecolor='red',ec='k', shrink=0.009),
#            horizontalalignment='left', verticalalignment='top',
#             );
#
#plt.annotate('Triv. cutoff', xy=(19, 0.5e-10), fontsize=31,color=[0.0,0.,0.9],  xycoords='data',
#             xytext=(0.64, 0.48), textcoords='axes fraction',
#             arrowprops=dict(facecolor=[0.,0.2,1],ec='k', shrink=0.009),
#             horizontalalignment='left', verticalalignment='top',
#             );


#plt.legend([p1, p2], [r'$\rm \,\,Trivial\,\,$', r'$\rm Topological$'],fontsize=28);

xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );
#plt.text(0.32, 5e-2, '(b) Circular Driver', fontsize=34, color="k");
#plt.text(2, 5e-2, '      Circular Driver', fontsize=31, color="k");
plt.yscale('log')
plt.xlim([3,29])

plt.grid(True)

plt.xlabel(r'$\rm \omega/\omega_0$', fontsize=37 );

ax.set_facecolor((0.1, 0.1, 0.1) )
fig.patch.set_facecolor('xkcd:mint green')
#plt.text(-.225, 0.5, r'$\rm I_{HHG}\,\,(arb.\,u.)$ ',
#         fontsize=37,
#         horizontalalignment='left',
#         verticalalignment='center',
#         rotation=90,
#         clip_on=False,
#         transform=plt.gca().transAxes)


plt.tick_params(labelsize = 37 );

#yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
yticks0  = [1e-20,1e-16, 1e-8,  1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.0e-3]);

#fname='HHG_Spectrum_Trivial_Vs_NonTrivial__Right_Circular_Polarization.pdf'
fname='Figure4bNew.pdf'
plt.savefig(fname, dpi = 400);

plt.show();
