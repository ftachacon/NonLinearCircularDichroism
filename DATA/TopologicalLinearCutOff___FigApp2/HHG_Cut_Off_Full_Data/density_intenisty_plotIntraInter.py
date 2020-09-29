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


plt.rc('lines', linewidth=2.5, color='k')


set_of_params   = sys.argv;
iparam             =  set_of_params[1];  #mag. flux phase

#################################
n 	     = 0;
ntemp    = 7;
Nt 	     = 101647#133411;
Nthalf 	 = int(Nt/30.);

if iparam=='Trivial':
    f 	 = open("filestrivial.dat","r");
else:
    f    = open("filesnontrivial.dat","r");


myinfo   = [];
phase 	 = [];
data 	 = [];
Intra_hhgMAP0  = []
Intra_hhgMAP   =np.zeros( (ntemp,Nthalf) )#[];

Inter_hhgMAP0  = []
Inter_hhgMAP   =np.zeros( (ntemp,Nthalf) )#[];

#################################
for line in f:
    myinfo.append(line);
    data.append(  np.loadtxt( line.strip() ) );
    phase.append( data[n][0,1] );
    Intra_hhgMAP0.append( data[n][1:Nthalf+1,2] + data[n][1:Nthalf+1,4] );
    Inter_hhgMAP0.append( data[n][1:Nthalf+1,1] + data[n][1:Nthalf+1,3] );
    n+=1;
    print("\nfile = ",line," n = ",n);

for i in range(n):
    Inter_hhgMAP[i,:] = Inter_hhgMAP0[i][:]
    Intra_hhgMAP[i,:] = Intra_hhgMAP0[i][:]

om = data[0][1:Nthalf+1,0];


print( '\n')
print( len(data[0][:,0]), 'omega-max =  ', om[-1], ' len( om ) = ' ,len(om));
print( phase );
print( '\nNlist = ', n, '  Nthalf = ', Nthalf );



#################################
#################################
width   = 11;
hight   = width/1.62;
fig     = plt.figure( figsize=(width,hight) );
ax1     = fig.add_axes([0.135, 0.1, 0.887, 0.835]);
for axis1 in ['top','bottom','left','right']:
    ax1.spines[axis1].set_linewidth(3.0);

nmin    = 1e-18;
nmax    = 5.e-1;
levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );

colorMap0='seismic'###'RdBu_r'#'gnuplot2'#'CMRmap_r'#"'bwr'#'coolwarm'#'gist_rainbow'#    #'PRGn'#'bwr'#'gnuplot2'#'CMRmap_r'#'coolwarm'#

cmap = plt.get_cmap(colorMap0)


norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
e0    = np.array([0.002,0.003,0.004,0.0045,0.005,0.006,0.007]);

X,Y     = np.meshgrid( om, e0 );

print( '\nshape of X ',X.shape, ' shape of e0 = ', e0.shape)

fer = interp2d(om, e0, np.log10(Inter_hhgMAP), kind='cubic')
fra = interp2d(om, e0, np.log10(Intra_hhgMAP), kind='cubic')

xnew = np.arange( 0., max(om)+5.e-3, 5.e-3 )
ynew = np.arange( min(e0), max(e0)+.00025, .00025)

data_er =  fer(xnew,ynew)
data_ra =  fra(xnew,ynew)
Xn, Yn = np.meshgrid(xnew, ynew)


#pcolormesh contourf
cf = plt.pcolormesh( Xn, Yn, np.power(10.,data_ra)+.99e-18,
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );
#cf = plt.pcolormesh( X, Y, data1
#                    ,vmin=nmin
#                    ,vmax=nmax #,levels=levels
#                    ,cmap=cmap );
xticks0  = np.arange(1,150,4);
plt.xticks( xticks0 );


#plt.xlabel('Harmonic Order', fontsize=21);


cb = plt.colorbar( cf, ticks=np.power(10.,np.arange( -18, 1, 4. ) ) );#mappable
cb.ax.tick_params(labelsize=23)

if iparam=='Trivial':
    plt.text(28, 0.0022, "Trivial Intra", fontsize=30, color=[1,1,1])
    plt.text(37.0, 0.00658, "(c)", fontsize=30, color=[1,1,1])
    m0=(0.007-0.002)/(32.5-13)
    x=np.arange(1.,34,.01);
    y=m0*(x-13)+0.002
    plt.plot(x,y,color=[0.,0.8,0],lw=3,linestyle='-')
    plt.ylabel(r'$E_0\,\, {\rm (a.u.)}$', fontsize=27);
else:
    plt.text(28, 0.0022, "Topo. Intra", fontsize=30, color=[1,1,1])
    plt.text(37.0, 0.00658, "(d)", fontsize=30, color=[1,1,1])
    m0=(0.007-0.002)/(34-14)
    x=np.arange(1.,38,.01);
    y=m0*(x-14)+0.002
    plt.plot(x,y,color=[0.,0.8,0],lw=3,linestyle='-')
print( 'm0 = ', m0)
#cb.set_tick_params(labelsize=16)
#plt.xlim([0, 27])
#fig.colorbar(cax, ticks=np.arange( 1.e-16,1.e0,1.e2 ))


plt.tick_params(labelsize=23);
plt.ylim([min(e0),max(e0)])
plt.xlim([0, 41])

if iparam=='Trivial':
    fname = 'FigApp2c.png'
else:
    fname = 'FigApp2d.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig(fileNamePicture, dpi = 150);


####################################
## Inter --Band contribution  ##
########################
fig     = plt.figure( figsize=(width,hight) );
ax1     = fig.add_axes([0.135, 0.126, 0.887, 0.835])
for axis1 in ['top','bottom','left','right']:
    ax1.spines[axis1].set_linewidth(3.0);

nmin    = 1e-18;
nmax    = 5.e-1;
levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0='seismic'
cmap = plt.get_cmap(colorMap0)
norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );


cf = plt.pcolormesh( Xn, Yn, np.power(10.,data_er)+.99e-18,
                    norm  = LogNorm( vmin=nmin, vmax=nmax ),
                    cmap = cmap );
xticks0  = np.arange(1,150,4);
plt.xticks( xticks0 );


cb = plt.colorbar( cf, ticks=np.power(10.,np.arange( -18, 1, 4. ) ) );#mappable
cb.ax.tick_params(labelsize=23)

if iparam=='Trivial':
    plt.text(28, 0.0022, "Trivial Inter", fontsize=30, color=[1,1,1])
    plt.text(37.0, 0.00658, "(e)", fontsize=30, color=[1,1,1])
    m0=(0.007-0.002)/(32.5-13)
    x=np.arange(1.,34,.01);
    y=m0*(x-13)+0.002
    plt.plot(x,y,color=[0.,0.8,0],lw=3,linestyle='-')
    plt.ylabel(r'$E_0\,\, {\rm (a.u.)}$', fontsize=27);
    plt.xlabel(r'${\rm Harmonic\,\,Order}$', fontsize=27);
else:
    plt.text(28., 0.0022, "Topo. Inter", fontsize=30, color=[1,1,1])
    plt.text(37.0, 0.00658, "(f)", fontsize=30, color=[1,1,1])
    m0=(0.007-0.002)/(34-14)
    x=np.arange(1.,38,.01);
    y=m0*(x-14)+0.002
    plt.plot(x,y,color=[0.,0.8,0],lw=3,linestyle='-')
    plt.xlabel(r'${\rm Harmonic\,\,Order}$', fontsize=29);
print ('m0 = ', m0)

plt.tick_params(labelsize=23);
plt.ylim([min(e0),max(e0)])
plt.xlim([0, 41])

if iparam=='Trivial':
    fname = 'FigApp2e.png'
else:
    fname = 'FigApp2f.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig(fileNamePicture, dpi = 150);





plt.show()

