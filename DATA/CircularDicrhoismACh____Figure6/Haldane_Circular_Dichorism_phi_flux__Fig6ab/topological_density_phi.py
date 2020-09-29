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
n 	     = 0;
ntemp    = 34;
Nt 	     = 133411;
Nthalf 	 = int(133411/16./2.);

if iparam=="LCP":
    f 	     = open("files1.dat","r");
else:
    f        = open("files2.dat","r");


phase 	 = np.loadtxt('DATA__PHASE__AXIS__'     + iparam + '.dat');
om       = np.loadtxt('DATA__FREQ__AXIS__'      + iparam + '.dat');
hhgMAP0  = np.loadtxt('DATA__HHG__PHASE__MAP__' + iparam + '.dat');

n        = len(phase);
Nthalf   = len(om);
hhgMAP   = np.reshape(hhgMAP0[:,2],(n,Nthalf));

print( phase );
print( '\nNlist = ', n, '  Nthalf = ', Nthalf );

HHG_MAP = hhgMAP #np.reshape( hhgMAP, (n,Nt-1) );



#################################
width   = 12;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
if iparam=='RCP':
    ax0=plt.axes([0.07, 0.125, 0.79, 0.825])
else:
    ax0=plt.axes([0.11, 0.125, 0.91, 0.825])
#ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
#nmin    = np.log10(1e-16);
#nmax    = np.log10(1e-2);

for axis in ['top','bottom','left','right']:
    ax0.spines[axis].set_linewidth(3.5);

nmin    = 1e-16;
nmax    = 1e-0;

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
colorMap0='CMRmap'#'RdBu_r'#'gnuplot2'#'bwr'#'coolwarm'#'gist_rainbow'#    #'seismic'#'gray_r'#'RdGy'#'PRGn'#'bwr'#'gnuplot2'#'CMRmap_r'#'coolwarm'#
cmap = plt.get_cmap(colorMap0);

norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
X,Y     = np.meshgrid(om,phase);

print ('shape of X ',X.shape);

f       = interp2d(om, phase, np.log10(HHG_MAP+1e-17), kind='cubic')
xnew = np.arange( 0., max(om), float(om[1]-om[0]) )
ynew = np.arange( 0., max(phase), .01)

data1 =  f(xnew,ynew)
#data1=np.log10(HHG_MAP+1e-17)
Xn, Yn = np.meshgrid(xnew, ynew)


cf = plt.pcolormesh( Xn, Yn, np.power(10.,data1),
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );
xticks0  = np.arange(1,40,4);
if iparam=='RCP':
    xticks0  = np.arange(1,32,4);

plt.xticks( xticks0 );
plt.xlabel(r'$\rm Harmonic\,\, Order$', fontsize=31);

if iparam=='LCP':
    cb = plt.colorbar( cf, ticks=[1.e-16,1.e-12,1.e-8,1.e-4,1.e-0] );
    cb.ax.tick_params(labelsize=25);

plt.tick_params(labelsize=25);
plt.ylim([0,3.12]);
plt.xlim([0, 33.]);

if iparam=='RCP':
    plt.text(26.5, 2.92, '(d)  ' + iparam, fontsize=30, color=[1,1,1]);
else:
    plt.text(26.5, 2.92, '(c)  ' + iparam, fontsize=30, color=[1,1,1]);

for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(2.5);

if iparam=='LCP':
    plt.ylabel(r'$\phi_0\, {\rm (rad.)}$', fontsize=31);
else:
    plt.ylabel("", fontsize=31);

if iparam=='RCP':
    chernNo = np.loadtxt('Phi__Vs__Chern__Num.dat');
    gapPhi = np.loadtxt('Phi__Vs__Energy__Gap.dat');
    
    ax1=plt.axes([0.875, 0.125, 0.112, 0.825]);
    ax1.plot(gapPhi[:,1],gapPhi[:,0],lw=4);
    ax1.set_yticks([]);
    ax1.set_xticks([0,0.1]);
    ax1.set_xlim([0.0,max(gapPhi[:,1])*1.04]);
    ax1.set_ylim([0,3.12]);
    ax1.tick_params(labelsize=20);
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2.5);
    shift=.0
    x1  = 0.52+shift;
    x2  = 2.630-shift;
    Npoints=3000
    y   = np.linspace(0.,np.pi,Npoints,endpoint=True)
    Npoints=len(y)
    cn  = [];
    cn0 = [];
    h=0.13;
    for i in range(Npoints):
        if (y[i]>=x1 and y[i]<=x2):
            cn.append(h)
        #cn0.append(0)
        else:
            cn.append(0)
            
        cn0.append(0.13)
    
    ax1.fill_between( cn,  0, y, color=[0.7,0.2,0.], alpha=0.150 )
    #ax1.fill_between( cn0, y/max(y)*x1,0, color=[0.2,0.2,0.2], alpha=0.950 )
    ax1.plot([0,0.13],[x1,x1],color=[0.,0.75,0],lw=8)
    ax1.plot([0,0.13],[x2,x2],color=[0.,0.75,0],lw=8)
    plt.text(0.0037, 2.92, "(e)", fontsize=30, color="k");
    plt.text(0.0037, 1.51, r'$C=+1$', fontsize=21, color="k");
    plt.text(0.042, 2.69,  r'$C=0$', fontsize=21, color="k");
    plt.text(0.042, 0.35,  r'$C=0$', fontsize=21, color="k");
    ax1.set_xlabel(r'$ E_{g}\,(a.u.)$ ', fontsize=25)

#print cn0



fname               = 'HHG__AND__MAP__'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig( fileNamePicture, dpi = 250 );




plt.show()
