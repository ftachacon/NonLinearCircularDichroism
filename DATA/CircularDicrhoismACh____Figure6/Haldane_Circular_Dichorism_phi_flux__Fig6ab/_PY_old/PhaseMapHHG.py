#!/usr/bin/python
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import os
import sys
import shutil
import cmath
import errno
from scipy.fftpack import fft, ifft

f      = open("files1.dat","r");
myinfo   = [];
phase      = [];
data      = [];
hhgMAP      = [];
n      = 0;
Nt      = 133412;
for line in f:
    myinfo.append(line);
    data.append(  np.loadtxt( line.strip() ) );
    phase.append( data[n][0,1] );
    hhgMAP.append( data[n][1:Nt,5] );
    print("\nfile = ",line);
    n+=1;


print( len(data[0][:,0]) );
print( phase );
print( '\nNlist = ', n, 'End of python test\n' );

HHG_MAP = hhgMAP #np.reshape( hhgMAP, (n,Nt-1) );

print (HHG_MAP.shape)



X,Y     = np.meshgrid(phase,kx);

nmin    = fparams[4]#occupationCB_nc.min();
nmax    = fparams[5]#occupationCB_nc.max();


levels = MaxNLocator(nbins=18).tick_values( nmin, nmax );


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.

colorMap0='bwr'#'coolwarm'#'gist_rainbow'#'RdBu'
    #'seismic'#'gray_r'#'RdGy'#'PRGn'#'bwr'#'gnuplot'#'gnuplot2'#'CMRmap_r'#'coolwarm'#'CMRmap_r'
cmap = plt.get_cmap(colorMap0)


norm = BoundaryNorm(levels
                    ,ncolors=cmap.N
                    ,clip=True
                    );


#####################################################
## Visualization
fontz  = 30;

width = 12
hight = width/1.62

fig = plt.figure(figsize=(width,hight));
ax1 = fig.add_axes([ax1, ax2, 0.85, 0.85])


#cf = plt.contourf( T
#                  ,Q
#                  ,occupationCB_nc
#                  ,levels=levels
#                  ,cmap=cmap );
#
#pcolormesh
if scale=="linear":
        cf=plt.pcolormesh( X
                      ,Y
                      ,occupationCB_nc
                      ,vmin=nmin
                      ,vmax=nmax #,levels=levels
                      ,cmap=cmap );
        
                      
if scale=="log":
            cf=plt.pcolormesh( X
                                ,Y
                                ,occupationCB_nc
                                ,norm = LogNorm( vmin=nmin
                                ,vmax=nmax )
                                ,cmap=cmap );
 
 
cb = plt.colorbar(cf);#mappable
plt.title(atitle,fontsize=21 )

plt.xlabel(xlab, fontsize=21);
plt.ylabel(ylab, fontsize=21);

plt.tick_params(labelsize=21);
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
