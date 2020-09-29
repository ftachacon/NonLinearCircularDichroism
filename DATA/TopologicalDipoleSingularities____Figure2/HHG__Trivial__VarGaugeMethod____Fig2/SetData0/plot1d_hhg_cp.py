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

ndat0='HHG__Phi0__0.060__M0t2__2.540__e0__1.000__G1.dat'
ndat1='HHG__Phi0__0.060__M0t2__2.540__e0__1.000__G2.dat'
ndat2='HHG__Phi0__0.060__M0t2__2.540__e0__1.000__G1G2.dat'

htrivial = np.loadtxt(ndat0);
htopo    = np.loadtxt(ndat1);
htri     = np.loadtxt(ndat2);

Ntime           = int(htrivial[0,0]);
dt              =   htrivial[0,5];

print "Ntime = ", Ntime, ";  dt = ",dt;


om          = htrivial[1:Ntime,0];
hhgtriv     = htrivial[1:Ntime,5];
hhgtopo     = htopo[1:Ntime,5];
htri1       = htri[1:Ntime,5]

#################################
#################################
width   = 11;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2.0);

hzero   = 0.25e-16;

p1,=plt.plot( om, hhgtriv + hzero, 'b',  lw=3.2    );
p2,=plt.plot( om, hhgtopo + hzero, '-r', lw=1.5 );
p3,=plt.plot( htri[1:Ntime,0], htri1 + hzero, '--g', lw=1.2 );

#plt.fill_between( om, hhgtopo, hhgtriv, color=[0.4,0.,0.95], alpha=0.250 )
#h1 = 0.50e-6
#h2 = 0.25e-10
#plt.plot([0,16.25],[h1,h1],color=[0.85,0.,0],lw=4,linestyle='--')
#plt.plot([0,19.2],[h2,h2],color=[0.,0.,0.95],lw=4,linestyle='--')

plt.legend([p1, p2,p3], [r'$\rm Gauge\,\,A\,\,$', r'$\rm Gauge\,\,B$',r'$\rm Variable-Gauge\,A-B$'],fontsize=28);

xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );
plt.text(1.25, 5e-2, "(a)", fontsize=34, color="k");
#plt.text(2, 5e-2, '      Circular Driver', fontsize=28, color="k");
plt.yscale('log')
plt.xlim([0,33])

plt.grid(True)

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=34 );
plt.ylabel(r'$\rm I_{HHG}(arb. u.)$', fontsize=34 );
plt.tick_params(labelsize = 30 );

yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([hzero*0.9,1.0e0]);

fname='Gauge__Comparison__Trivial.pdf'
plt.savefig(fname, dpi = 400);

plt.show();
