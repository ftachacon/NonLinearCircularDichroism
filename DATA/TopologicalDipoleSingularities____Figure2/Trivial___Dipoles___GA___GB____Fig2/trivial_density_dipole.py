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

#Input parameter: this is provided by terminal
set_of_params   = sys.argv;
iparam          =  set_of_params[1];  #mag. flux phase

#################################
#
kx 	 = np.loadtxt('kxAxis.dat');
ky   = np.loadtxt('kyAxis.dat');

temp = 'RealDipoleTrivialGauge'+iparam+'.dat';

Dip1 = np.transpose( np.loadtxt( temp ) );
print ('\nSize of Dip1 = ', np.shape( Dip1 ));

X,Y     = np.meshgrid(kx,ky);

nmin    = -2#Dip1.min();
nmax    = 2.01#Dip1.max();

#################################
width   = 10.2;
hight   = width#/1.62;
fig     = plt.figure( figsize=(width,hight) );
#           x   y      width  high
AxConf1 = [ 0.106, 0.116, 0.88,   0.87]
AxConf2 = [ 0.10, 0.116, 0.887,   0.87]
ax0     = plt.axes(AxConf2)#[0.11, 0.115, 0.865, 0.87]);
for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(3);

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0='seismic'#'CMRmap'#'RdBu_r'#'gnuplot2'#'bwr'#'coolwarm'#'gist_rainbow'#


cmap = plt.get_cmap(colorMap0);
norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );


f = interp2d(kx, ky, Dip1, kind='cubic');

xnew = np.arange( min(kx), max(kx), max(kx)*2/501 )
ynew = np.arange( min(ky), max(ky), max(ky)*2/501)


Xn, Yn = np.meshgrid(xnew, ynew)
Dat2    = f(xnew,ynew)

cf = plt.pcolormesh( Xn, Yn, Dat2,
                 vmin=nmin,
                 vmax=nmax,
                 cmap = cmap );

plt.tick_params(labelsize=40);
xticks0  = np.arange(-2,2,1);
plt.xticks( xticks0 );
plt.xlim([min(kx),max(kx)]);
plt.ylim([min(ky),max(ky)]);

K1 = np.array([-1.279333,0])
K2 = np.array([1.279333,0])
K3 = np.array([-1.279333/2,1.1079])
K4 = np.array([1.279333/2,1.1079])
K5 = np.array([-1.279333/2,-1.1079])
K6 = np.array([1.279333/2,-1.1079])

gren=[0,1,0]

plt.plot([K1[0],K3[0]],[K1[1],K3[1]], color=gren,lw=3)
plt.plot([K4[0],K2[0]],[K4[1],K2[1]], color=gren,lw=3)
plt.plot([K6[0],K2[0]],[K6[1],K2[1]], color=gren,lw=3)
plt.plot([K3[0],K4[0]],[K3[1],K4[1]], color=gren,lw=3)
plt.plot([K1[0],K5[0]],[K1[1],K5[1]], color=gren,lw=3)
plt.plot([K5[0],K6[0]],[K5[1],K6[1]], color=gren,lw=3)

wid=21
plt.plot(K1[0],K1[1],'wo',lw=6,markersize=wid)
plt.plot(K6[0],K6[1],'wo',lw=6,markersize=wid)
plt.plot(K4[0],K4[1],'wo',lw=6,markersize=wid)

plt.plot(K5[0],K5[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K3[0],K3[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K2[0],K2[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)

xticks0  = np.arange( -2. ,2.2, 1);
yticks0  = np.arange( -2. ,2.2, 1);
plt.xticks( xticks0 );
plt.yticks( yticks0 );
plt.xlim([min(kx)*1.045,max(kx)*1.045]);
plt.ylim([min(ky),max(ky)]);
print ("kx = ", min(kx))

#ax1.set_xticks( xticks0 );


plt.xlabel(r'$\rm k_x\,\, (\,a.u.)$', fontsize=37);
if iparam=='1':
    #plt.ylabel(r'$\rm k_y\,\, (a.u.)$', fontsize=32,multialignment='center');
    plt.text(-2.48,0.3,r'$\rm k_y\,\, (a.u.)$', fontsize=37,rotation=90);
    plt.text(-1.94, 1.97, r'$\rm (a) \,\, Gauge\, A$'
             , fontsize=37, color='k'
             , bbox=dict(boxstyle="round"
             , ec=(1, 1.0, 1.0), fc=(1, 1, 1)) )

    plt.text( -1.1, -.075, r"$\rm K'$", fontsize=27, color='k',bbox=dict(boxstyle="round"
                , ec=(0.0, 1.0, 1.0), fc=(1, 1, 1)) );
    
    #plt.text( -1.0, -.1, r"$\rm K'$", fontsize=27, color='k',bbox=dict(boxstyle="round"
    #    , ec=(0.0, 1.0, 1.0), fc=(1, 1, 1)) );
    
    plt.text( -0.63, -.86, r"$\rm K$", fontsize=27, color='k',bbox=dict(boxstyle="round"
        ,ec=(1.0, 1.0, 0.0), fc=(1, 1, 1)) );


#plt.text(4.5, 0.5, 'text -45', props, rotation=-45)

             #ax0.yaxis.set_label_position("left")
#ax0.yaxis.tick_left()
#plt.text(-1.82, 1.50, r"$\rm  Trivial\,\, Dipole $", fontsize=26, color='k')

if iparam=='2':
    #plt.ylabel(r'$\rm k_y\, (a.u.)$', fontsize=33);
    plt.text(-1.94, 1.97, r'$\rm (b) \,\, Gauge\, B$', fontsize=37, color='k'
        ,bbox=dict(boxstyle="round", ec=(1, 1.0, 1.0), fc=(1, 1, 1)) )
#plt.text(-1.82, 1.50, r"$\rm  Trivial\,\, Dipole $", fontsize=26, color='k')

fname               = 'DipoleGauge'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig( fileNamePicture, dpi = 150 );

plt.show()
