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
iparam          =  set_of_params[1];  #mag. flux phase

#################################
#
kx 	 = np.loadtxt('kxAxis.dat');
ky   = np.loadtxt('kyAxis.dat');

Dip1 = np.transpose(np.loadtxt('RealDipoleTopoGauge'+iparam+'.dat'));
#Dip2 = np.loadtxt('RealDipoleTopoGauge2.dat');
#Dip3 = np.loadtxt('RealDipoleTopoGauge3.dat');

X,Y     = np.meshgrid(kx,ky);

nmin    = -2.0#Dip1.min();
nmax    = 2.01#Dip1.max();

#################################
width   = 11;
hight   = width/1.62;
fig     = plt.figure( figsize=(width,hight) );

AxConf1=[0.078, 0.115, 0.91, 0.87]
AxConf2=[0.065, 0.115, 0.92, 0.87]
ax0     = plt.axes(AxConf1);

for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(3);

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0='seismic'#'hot'#'RdBu_r'#'CMRmap'#'gnuplot2'#'bwr'#''gist_rainbow'#
cmap = plt.get_cmap(colorMap0);
norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );

f = interp2d( kx, ky, Dip1, kind='cubic' );

xnew = np.arange( min(kx), max(kx), max(kx)*2/501 )
ynew = np.arange( min(ky), max(ky), max(ky)*2/501)

Xn, Yn = np.meshgrid(xnew, ynew)
Dat2    = f(xnew,ynew)

cf = plt.pcolormesh( Xn, Yn, Dat2,
                 vmin=nmin,
                 vmax=nmax,
                 cmap = cmap );

plt.tick_params(labelsize=40);

K1 = np.array([-1.279333,0]);
K2 = np.array([1.279333,0]);
K3 = np.array([-1.279333/2,1.1079]);
K4 = np.array([1.279333/2,1.1079]);
K5 = np.array([-1.279333/2,-1.1079]);
K6 = np.array([1.279333/2,-1.1079]);

gren=[0,1,0]

plt.plot([K1[0],K3[0]],[K1[1],K3[1]],color=gren,lw=3)
plt.plot([K4[0],K2[0]],[K4[1],K2[1]],color=gren,lw=3)
plt.plot([K6[0],K2[0]],[K6[1],K2[1]],color=gren,lw=3)
plt.plot([K3[0],K4[0]],[K3[1],K4[1]],color=gren,lw=3)
plt.plot([K1[0],K5[0]],[K1[1],K5[1]],color=gren,lw=3)
plt.plot([K5[0],K6[0]],[K5[1],K6[1]],color=gren,lw=3)

wid=21
plt.plot(K1[0],K1[1],'wo',lw=6,markersize=wid)
plt.plot(K6[0],K6[1],'wo',lw=6,markersize=wid)
plt.plot(K4[0],K4[1],'wo',lw=6,markersize=wid)
plt.plot(K5[0],K5[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K3[0],K3[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K2[0],K2[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)

#cb = plt.colorbar( cf, ticks=np.power(10.,np.arange( -18, 1, 4. ) ) );#mappable
#cb = plt.colorbar( cf, ticks=np.arange(-2,2,1) );
#cb.ax.tick_params(labelsize=23)
xticks0  = np.arange( -2. ,2.2, 1);
yticks0  = np.arange( -2. ,2.2, 1);
plt.xticks( xticks0 );
plt.yticks( yticks0 );

plt.xlim([min(kx)*1.045,max(kx)*1.045]);
plt.ylim([min(ky),max(ky)]);
print ("kx = ", min(kx))


plt.xlabel(r'$\rm k_x\,\, (\,a.u.)$', fontsize=37);
if iparam=='1':
    #plt.ylabel(r'$\rm k_y\, (a.u.)$', fontsize=33);
    plt.text( -1.94, 1.97, r'$\rm (c) \,\, Gauge\, A$', fontsize=37, color='k',bbox=dict(boxstyle="round", ec=(1, 1.0, 1.0), fc=(1, 1, 1)) );
    
    plt.text( -1.1, -.075, r"$\rm K'$", fontsize=27, color='k',bbox=dict(boxstyle="round"
       , ec=(0.0, 1.0, 1.0), fc=(1, 1, 1)) );
        
        #plt.text( -1.0, -.1, r"$\rm K'$", fontsize=27, color='k',bbox=dict(boxstyle="round"
        #    , ec=(0.0, 1.0, 1.0), fc=(1, 1, 1)) );

    plt.text( -0.63, -.86, r"$\rm K$", fontsize=27, color='k',bbox=dict(boxstyle="round"
        , ec=(1.0, 1.0, 0.0), fc=(1, 1, 1)) );

#plt.text( -1.82, 1.50, r"$\rm  Topological\,\, Dipole $", fontsize=26, color='k')

if iparam=='2':
    #plt.ylabel(r'$\rm k_y\, (a.u.)$', fontsize=33);
    plt.text(-1.94, 1.97, r'$\rm (d) \,\, Gauge\, B$', fontsize=37,
             color='k',
             bbox=dict(boxstyle="round"
            ,ec=(1, 1.0, 1.0), fc=(1, 1, 1)) );


#plt.text(-1.82, 1.50, r"$\rm  Topological\,\, Dipole $", fontsize=26, color='k')

#plt.tick_params(labelsize=25);
#if iparam=='RCP':
#    plt.text(26.5, 2.92, '(d)  ' + iparam, fontsize=30, color=[1,1,1]);
#else:
#    plt.text(26.5, 2.92, '(c)  ' + iparam, fontsize=30, color=[1,1,1]);


if iparam == '1':
    ldip1 = np.loadtxt('RealDipoleTopoLineskxGauge'+iparam+'.dat');
    ax1=plt.axes([0.12, 0.174, 0.34, 0.134]);
    ax1.plot(ldip1[:,0],ldip1[:,1],linestyle='-',lw=1.25,color=[0,0,0],label='kx = -1.85');
    ax1.plot(ldip1[:,0],ldip1[:,2],lw=1.25,color=[0,1,0],label='kx = -1.80');
    ax1.plot(ldip1[:,0],ldip1[:,3],lw=1.25,color=[00,0,1],label='kx = -1.75');
    ax1.plot(ldip1[:,0],ldip1[:,4],lw=1.8,color=[1,0,0],label='kx = -1.70');
    #ax1.tick_params(axis='x', colors='y')

#ax1.legend(bbox_to_anchor=(0.1, 0.1))
    ax1.grid(True)
    ax1.set_xticks([-2,-1,0,1,2]);
    ax1.set_yticks([-2,0,2]);
    ax1.set_xlim([min(ky),max(ky)]);
    ax1.set_ylim([-2,2]);
    ax1.set_xlabel(r'$\rm \, ky\,\,(a.u.) \,\,\,\,\,\,\,\,$',color='g', fontsize=14);

    ax1.set_xticks( xticks0 );
    gren=[0,1,0]
    ax1.spines['bottom'].set_color(gren)
    ax1.spines['top'].set_color(gren)
    ax1.spines['left'].set_color(gren)
    ax1.spines['right'].set_color(gren)
    ax1.xaxis.label.set_color(gren)
    ax1.tick_params(axis='x', colors=gren)
    ax1.tick_params(axis='y', colors=gren)

#   ax1.yaxis.label.set_color([0,1,0.2])
#   ax1.xaxis.label.set_color([0,1,0.2])
    ax1.tick_params(labelsize=14,color=[0.,1,0.2]);

#for axis in ['top','bottom','left','right']:
#ax1.spines[axis].set_linewidth(2.5);

if iparam == '2':
    ldip1 = np.loadtxt('RealDipoleTopoLineskyGauge'+iparam+'.dat');
    #ax1=plt.axes([0.14, 0.165, 0.34, 0.14]);
    ax1=plt.axes([0.12, 0.174, 0.34, 0.134]);
    ax1.plot(ldip1[:,0],ldip1[:,1],linestyle='-',lw=1.25,color=[0,0,0],label='kx = -1.85');
    ax1.plot(ldip1[:,0],ldip1[:,2],lw=1.25,color=[0,1,0],label='kx = -1.80');
    ax1.plot(ldip1[:,0],ldip1[:,3],lw=1.25,color=[00,0,1],label='kx = -1.75');
    ax1.plot(ldip1[:,0],ldip1[:,4],lw=1.8,color=[1,0,0],label='kx = -1.70');
    
    #ax1.legend(bbox_to_anchor=(0.1, 0.1))
    ax1.grid(True)
    ax1.set_xticks([-2,-1,0,1,2]);
    ax1.set_yticks([-2,0,2]);
    ax1.set_xlim([min(kx),max(kx)]);
    ax1.set_ylim([-2,2]);
    ax1.set_xlabel(r'$\rm \, kx\,\,(a.u.) \,\,\,\,\,\,\,\,$',color='g', fontsize=14);
#.tick_params(labelsize=12);
    ax1.tick_params(labelsize=16,color=[0.,1,0.2]);

    ax1.spines['bottom'].set_color(gren)
    ax1.spines['top'].set_color(gren)
    ax1.spines['left'].set_color(gren)
    ax1.spines['right'].set_color(gren)
    ax1.xaxis.label.set_color(gren)
    ax1.tick_params(axis='x', colors=gren)
    ax1.tick_params(axis='y', colors=gren)


fname               = 'DipoleGauge'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig( fileNamePicture, dpi = 150 );

plt.show()
