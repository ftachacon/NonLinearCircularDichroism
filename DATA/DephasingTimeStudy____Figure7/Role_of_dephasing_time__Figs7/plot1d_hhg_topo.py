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

temp0='HHG__Phi0__1.160__M0t2__2.540__e0__-1.000__T2__220.000__.dat'
temp1='HHG__Phi0__1.160__M0t2__2.540__e0__-1.000__T2__330.000__.dat'
temp2='HHG__Phi0__1.160__M0t2__2.540__e0__-1.000__T2__440.000__.dat'

htrivial0 = np.loadtxt(temp0);
htrivial1 = np.loadtxt(temp1);
htrivial2 = np.loadtxt(temp2);

T2          = np.array([220,330,440])

Ntime           = int(htrivial0[0,0]);
dt              = htrivial0[0,5];

print( "Ntime = ", Ntime, ";  dt = ",dt);

om           = htrivial0[1:Ntime,0];
hhgtriv0     = htrivial0[1:Ntime,5];
hhgtriv1     = htrivial1[1:Ntime,5];
hhgtriv2     = htrivial2[1:Ntime,5];

#################################
#################################
width   = 11;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
#ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
#for axis in ['top','bottom','left','right']:
#    ax.spines[axis].set_linewidth(3.0);

ax      = fig.add_axes([0.165, 0.14, 0.827, 0.84]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);

hzero   = 0.25e-16;

p1,=plt.plot( om, hhgtriv0 + hzero, 'b', lw=2   );
p2,=plt.plot( om, hhgtriv1 + hzero, 'g', lw=1.5   );
p3,=plt.plot( om, hhgtriv2 + hzero, 'r', lw=1.0   );

#plt.fill_between( om, hhgtopo, hhgtriv, color=[0.4,0.,0.95], alpha=0.250 )
#plt.legend([p1, p2], [ r'$C\,=\,\,\,\,0$', r'$C\,=\,+1$'],fontsize=28);

leg0=r'$T_2 \,=\, $' +str('%.0f'%T2[0]) + r'  $\rm a.u.$'
leg1=r'$\,\,\,\,=\, $' +str('%.0f'%T2[1])
leg2=r'$\,\,\,\,=\, $' +str('%.0f'%T2[2])

plt.legend([p1, p2,p3], [leg0 , leg1,leg2],fontsize=27);
xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

plt.text(2, 5e-1, r'$\rm(b)\,RCP\,\,Topological $', fontsize=31, color="k");


plt.yscale('log')
plt.xlim([0,35])

plt.grid(True)



#plt.xlabel(r'$\rm Harmonic-Order$', fontsize=34 );


yticks0  = [1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4 ]#10., 15, 5.0];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.0e1])
plt.ylabel(r'$\rm I_{HHG}\,\, (\,arb.\, u.)$', fontsize=32 );
plt.tick_params(labelsize = 31 );

#fname='Dephasing__Dependence__HHG__RCP__Topological.pdf'
fname='Figure7b.pdf'
plt.savefig(fname, dpi = 400);

plt.show()



