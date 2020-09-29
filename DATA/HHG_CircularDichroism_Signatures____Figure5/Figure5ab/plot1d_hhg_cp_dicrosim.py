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

set_of_params   = sys.argv;
iparam         =  set_of_params[1];  #mag. flux phase param

plt.rc('lines', lw=3, color='k')


if iparam=='Trivial':
    rcp = np.loadtxt("HHG__Phi0__0.060__M0t2__2.540__e0__-1.000__.dat");
    lcp = np.loadtxt("HHG__Phi0__0.060__M0t2__2.540__e0__1.000__.dat");
else:
    rcp = np.loadtxt("HHG__Phi0__1.160__M0t2__2.540__e0__-1.000__.dat");
    lcp = np.loadtxt("HHG__Phi0__1.160__M0t2__2.540__e0__1.000__.dat");

Ntime           = int(rcp[0,0]);
dt              =   rcp[0,5];

#print "Ntime = ", Ntime, ";  dt = ",dt;


om          = rcp[1:Ntime,0];
hhg__rcp     = rcp[1:Ntime,5];
hhg__lcp     = lcp[1:Ntime,5];


#################################
#################################
width   = 11.4;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);

hzero   = 0.25e-16;

p1,=plt.plot( om, hhg__rcp + hzero, 'b', lw=2   );
p2,=plt.plot( om, hhg__lcp + hzero, '-r', lw=1.1 );

plt.fill_between( om, hhg__lcp, hhg__rcp, color=[0.24,0.9,0.0], alpha=0.250 )


#plt.legend([p1, p2], [ r'$C\,=\,\,\,\,0$', r'$C\,=\,+1$'],fontsize=28);

plt.legend([p1, p2], [r'$\rm RCP\,\,\,\,$', r'$\rm LCP $'],fontsize=37);

xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

if iparam=='Trivial':
    plt.text(0.32, 5e-2, '(a) Triv. Phase', fontsize=37, color="k");
#plt.text(2, 5e-2, '    Trivial Phase', fontsize=28, color="k");
else:
    plt.text(0.32, 5e-2, '(b) Topo. Phase', fontsize=37, color="k");
#plt.text(2, 5e-2, '    Topological Phase', fontsize=28, color="k");

plt.yscale('log')
plt.xlim([0,33])

plt.grid(True)

plt.xlabel(r'$\rm Harmonic\,Order$', fontsize=34 );

if iparam=='Trivial':
    plt.ylabel(r'$\rm I_{HHG}$ (arb. u.)', fontsize=34 );


plt.tick_params(labelsize = 30 );

yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.0e0])

if iparam=='Trivial':
    fname='LCP_Vs_RCP_Trivial.pdf'
else:
    fname='LCP_Vs_RCP_NonTrivial.pdf'

plt.savefig(fname, dpi = 400);

plt.show()



