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
iparam             =  set_of_params[1];  #mag. flux phase

#################################
Mt2 	 = np.loadtxt('DATA__mt2__AXIS__'     + iparam + '.dat');
om       = np.loadtxt('DATA__FREQ__AXIS__'    + iparam + '.dat');
hhgMAP0  = np.loadtxt('DATA__HHG__mt2__MAP__' + iparam + '.dat');


#################################
n        = len(Mt2);
Nthalf   = len(om);
hhgMAP   = np.reshape(hhgMAP0[:,2],(n,Nthalf));

print( Mt2 );
print( '\nNlist = ', n, '  Nthalf = ', Nthalf );


HHG_MAP = hhgMAP #np.reshape( hhgMAP, (n,Nt-1) );

#################################
width   = 12;
hight   = width/1.62;
fig     = plt.figure( figsize=(width,hight) );
ashift=-0.03
if iparam=='RCP':
    ax0=plt.axes([0.10+ashift, 0.13, 0.73, 0.835]);
else:
    ax0=plt.axes([0.11, 0.13, 0.91, 0.835]);


for axis in ['top','bottom','left','right']:
    ax0.spines[axis].set_linewidth(3.5);


hzero   = 1.e-20
nmin    = hzero;
nmax    = 1e-0;

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0='CMRmap'
cmap = plt.get_cmap(colorMap0);
norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
X,Y     = np.meshgrid(om,Mt2);

print( 'shape of X =', X.shape);

f       = interp2d(om, Mt2, np.log10(HHG_MAP+hzero), kind='cubic')
xnew = np.arange( 0., max(om), float(om[1]-om[0]) )
ynew = np.arange( 0., max(Mt2), .01)

data1   =  f(xnew,ynew);
Xn, Yn  = np.meshgrid(xnew, ynew);


cf = plt.pcolormesh( Xn, Yn, np.power(10.,data1),
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );
xticks0  = np.arange(1,40,4);
if iparam=='RCP':
    xticks0  = np.arange(1,35,4);

plt.xticks( xticks0 );
plt.xlabel(r'$\rm Harmonic\,\, Order$', fontsize=31);
yaxis = np.arange( -20, 0, 4. );

if iparam=='LCP':
    cb = plt.colorbar( cf, ticks=np.power(10.,yaxis ) );
    cb.ax.tick_params(labelsize=25);

plt.tick_params(labelsize=25);
plt.ylim([0,max(Mt2)]);
plt.xlim([0, 33.]);

labelpos2=8.4
if iparam=='RCP':
    plt.text(25.5, labelpos2, '(e)  ' + iparam, fontsize=30, color=[1,1,1]);
else:
    plt.text(25.1, labelpos2, '(d)  ' + iparam, fontsize=30, color=[1,1,1]);

for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(2.5);

if iparam=='LCP':
    plt.ylabel(r'$ M_0/t_2$', fontsize=31);
else:
    plt.ylabel("", fontsize=31);

if iparam=='RCP':
    gapPhi  = np.loadtxt('HPD__Param.dat');
    eg1     = gapPhi[:,3]
    mt2     = gapPhi[:,1]
    ax1=plt.axes([0.875+ashift, 0.13, 0.112, 0.835]);
    ax1.plot(eg1,mt2,lw=2);
    ax1.set_yticks([]);
    ax1.set_xticks([0,0.1,0.2]);
    ax1.set_xlim([0.0,max(eg1)*1.05]);
    ax1.set_ylim([0,max(mt2)]);
    ax1.tick_params(labelsize=20);
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2.5);
    shift   = 0.0;
    x1      = 0.0;
    x2      = 3*np.sqrt(3)-shift;
    Npoints = 3000;
    y       = np.linspace(0.,max(mt2),Npoints,endpoint=True);
    Npoints = len(y);
    cn      = [];
    cn0     = [];
    h       = max(eg1);
    for i in range(Npoints):
        if (y[i]>=x1 and y[i]<=x2):
            cn.append(h);
        else:
            cn.append(0);
            
        cn0.append(h)
    
    ax1.fill_between( cn,  0, y, color=[0.0,0.8,0.99], alpha=0.04 )
    ax1.plot([0,h],[x2,x2],color=[0.,0.75,0],lw=8);
    plt.text(0.0037, labelpos2, "(f)", fontsize=30, color="k");
    plt.text(0.042, 4.4, r'$C=+1$', fontsize=20, color="k");
    plt.text(0.042, 5.7,  r'$C=\,0$', fontsize=21, color="k");
    ax1.set_xlabel(r'${ E}_g\,(a.u.)$ ', fontsize=27);

fname               = 'HHG__AND__MAP__'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig( fileNamePicture, dpi = 150 );

plt.show()
