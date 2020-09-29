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


#################################
n 	     = 0;
ntemp    = 7;
Nt 	     = 101647#133411;
Nthalf 	 = int(Nt/30.);

#f 	 = open("files1.dat","r");
f      = open("files2.dat","r");


myinfo   = [];
phase 	 = [];
data 	 = [];
hhgMAP0   = []
hhgMAP   =np.zeros( (ntemp,Nthalf) )#[];


#################################
for line in f:
    myinfo.append(line);
    data.append(  np.loadtxt( line.strip() ) );
    phase.append( data[n][0,1] );
    hhgMAP0.append( data[n][1:Nthalf+1,5] );
    n+=1;
    print("\nfile = ",line," n = ",n);

for i in range(n):
    hhgMAP[i,:] = hhgMAP0[i][:]

om = data[0][1:Nthalf+1,0];


print( '\n')
print( len(data[0][:,0]), 'omega-max =  ', om[-1], ' len( om ) = ' ,len(om));
print( phase );
print( '\nNlist = ', n, '  Nthalf = ', Nthalf );

HHG_MAP = hhgMAP #np.reshape( hhgMAP, (n,Nt-1) );



#################################
#################################
width   = 11;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
ax1     = fig.add_axes([0.125, 0.125, 0.89, 0.825])


#nmin    = np.log10(1e-16);
#nmax    = np.log10(1e-2);
nmin    = 1e-18;
nmax    = 1e-2;

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
colorMap0='gnuplot2'#'CMRmap_r'#"RdBu_r"#'bwr'#'coolwarm'#'gist_rainbow'#    #'seismic'#'gray_r'#'RdGy'#'PRGn'#'bwr'#'gnuplot2'#'CMRmap_r'#'coolwarm'#

cmap = plt.get_cmap(colorMap0)


norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
e0    = np.array([0.002,0.003,0.004,0.0045,0.005,0.006,0.007]);

X,Y     = np.meshgrid( om, e0 );

print( '\nshape of X ',X.shape, ' shape of e0 = ', e0.shape)

f = interp2d(om, e0, np.log10(HHG_MAP), kind='cubic')

xnew = np.arange( 0., max(om)+5.e-3, 5.e-3 )
ynew = np.arange( min(e0), max(e0)+.00025, .00025)

data1 =  f(xnew,ynew)
Xn, Yn = np.meshgrid(xnew, ynew)


#pcolormesh contourf
cf = plt.pcolormesh( Xn, Yn, np.power(10.,data1)+.99e-18,
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );

#cf = plt.pcolormesh( X, Y, data1
#                    ,vmin=nmin
#                    ,vmax=nmax #,levels=levels
#                    ,cmap=cmap );

xticks0  = np.arange(1,150,4);
plt.xticks( xticks0 );


plt.xlabel('Harmonic Order', fontsize=21);
plt.ylabel(r'$E_0 {\rm (a.u.)}$', fontsize=21);

cb = plt.colorbar( cf, ticks=np.power(10.,np.arange( -18, 1, 4. ) ) );#mappable

#cb = plt.colorbar( cf );#mappable
cb.ax.tick_params(labelsize=23)

plt.text(32.5, 0.0022, "Topological", fontsize=20, color=[1,1,1])
plt.text(37.0, 0.0065, "(b)", fontsize=25, color=[1,1,1])

#cb.set_tick_params(labelsize=16)
#plt.xlim([0, 27])
#fig.colorbar(cax, ticks=np.arange( 1.e-16,1.e0,1.e2 ))


plt.tick_params(labelsize=21);
plt.ylim([min(e0),max(e0)])
plt.xlim([0, 41])

fname = 'HHG__AND__MAP__LinearPol__Topological.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig(fileNamePicture, dpi = 150);


#############################
#ntotal = ntemp*Nthalf
#outdata=np.zeros( (ntotal,3))
#outname='DATA__HHG__PHASE__MAP__Trivial.dat'
outname='DATA__HHG__PHASE__MAP__Topological.dat'
afile   = open( outname, 'w' );

for i in range(ntemp):
    for j in range(Nthalf):
        temp   = str('%.4e'%e0[i]) + '     ' + str( '%.4e'%om[j] ) +'     ' + str( '%.16e'%HHG_MAP[i,j] );
        afile.write( temp );
        afile.write("\n")
afile.close()

plt.show()

