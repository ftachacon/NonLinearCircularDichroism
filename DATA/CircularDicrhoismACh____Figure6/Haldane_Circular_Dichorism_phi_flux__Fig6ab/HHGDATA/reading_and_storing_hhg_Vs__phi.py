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
from scipy.interpolate import interp1d

#################################
##
set_of_params   = sys.argv;
iparam          = set_of_params[1];

#################################
##
n 	     = 0;
ntemp    = 34;
Nt 	     = 133411;
Nthalf 	 = int(133411/16./2.);
print ('Nt" = ', Nthalf, '; Nt', Nt)


if iparam=="LCP":
    f 	     = open("files1.dat","r");
else:
    f        = open("files2.dat","r");

myinfo      = [];
phase 	    = [];
data 	    = [];
hhgMAP0     = [];


#################################
###
dwx1        = 0.025;
om1         = np.arange( dwx1, 40.+dwx1, dwx1 );
print( 'max(om1) = ',max(om1), ' min(om1) = ', min(om1));
NtNew1      = len(om1);
hzero       = 1.0e-16;

#################################
###
for line in f:
    myinfo.append(line);
    data.append(  np.loadtxt( line.strip() ) );
    phase.append( data[n][0,1] );
    spectra0    = data[n][1:Nthalf,5];
    xom0        = data[n][1:Nthalf,0];
    
    print( 'max(xom0) = ', max(xom0), ' min(xom0) = ',min(xom0) );
    
    f0          = interp1d( xom0, np.log10( spectra0 + hzero ), kind='slinear' );
    spectra1    = f0( om1 );
    hhgMAP0.append( spectra1 );
    n+=1;
    print( 'n = ', n, ';      file = ',line );

hhgMAP   = np.zeros( (n,NtNew1) );
for i in range(n):
    hhgMAP[i,:] = hhgMAP0[i][:];

print( len(data[0][:,0]), 'omega-max =  ', om1[-1] );
print( phase );
print( '\nNlist = ', n, '  Nthalf = ',NtNew1 );


#################################
nmin    = hzero; #np.log10(hzero);
nmax    = 1e-1;  #np.log10(1e-1);
width   = 11;
hight   = width/1.62;
fig     = plt.figure( figsize=(width,hight) );
plt.axes([0.07, 0.125, 0.79, 0.825]);
plt.rc('lines', lw=3, color='k');
levels      = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0   ='CMRmap'#'RdBu_r'#'gnuplot2'#'bwr'#'coolwarm'#'gist_rainbow'#
cmap        = plt.get_cmap(colorMap0);
norm        = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
X,Y         = np.meshgrid(om1,phase);
f           = interp2d( om1, phase, hhgMAP, kind='cubic' );
xnew        = om1;      #np.arange( min(om1), max(om1), float(om1[1]-om1[0]) );
ynew        = np.arange( 0., max(phase), .005);
data1       = f(xnew,ynew);
Xn, Yn      = np.meshgrid(xnew, ynew);

#pcolormesh contourf
#cf  = plt.pcolormesh( X, Y, np.power(10.,data1),
#                norm  = LogNorm( vmin=nmin, vmax=nmax ),
#                 cmap = cmap );
cf = plt.pcolormesh( Xn, Yn, np.power(10.,data1),
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );
#cf = plt.pcolormesh( Xn, Yn, data1
#                    ,vmin=nmin
#                    ,vmax=nmax #,levels=levels
#                    ,cmap=cmap );


xticks0  = np.arange(1,32,4);
plt.xticks( xticks0 );
plt.xlabel(r'$\rm Harmonic\,\, Order$', fontsize=31);

if iparam=='LCP':
    cb = plt.colorbar( cf );#mappable
    cb.ax.tick_params(labelsize=25);

plt.tick_params(labelsize=25);
plt.ylim([0,3.12]);
plt.xlim([0, 33.]);


if iparam=='LCP':
    plt.ylabel(r'$\phi_0 {\rm (rad.)}$', fontsize=31);
else:
    plt.ylabel("", fontsize=31);

if iparam=='RCP':
    chernNo = np.loadtxt('Phi__Vs__Chern__Num.dat');
    gapPhi  = np.loadtxt('Phi__Vs__Energy__Gap.dat');
    
    ax1=plt.axes([0.865, 0.125, 0.112, 0.825]);
    ax1.plot(gapPhi[:,1],gapPhi[:,0],lw=2);
    ax1.set_yticks([]);
    ax1.set_xticks([0,0.1]);
    ax1.set_xlim([0.0,max(gapPhi[:,1])]);
    ax1.set_ylim([0,3.12]);
    ax1.tick_params(labelsize=20);
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2.5);

fname               = 'HHG__AND__MAP__'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig(fileNamePicture, dpi = 150);
outname             = 'DATA__HHG__PHASE__MAP__' + iparam + '.dat'
afile               = open( outname, 'w' );

for i in range(n):
    for j in range(NtNew1):
        temp   = str('%.4e'%phase[i]) + '     ' + str( '%.4e'%om1[j] ) +'     ' + str( '%.16e'%np.power( 10.,hhgMAP[i,j]) );
        afile.write( temp );
        afile.write("\n");
    afile.write("\n");
afile.close();


outname='DATA__PHASE__AXIS__' + iparam + '.dat'
afile   = open( outname, 'w' );
for i in range(n):
    temp   = str('%.4e'%phase[i]) ;
    afile.write( temp );
    afile.write( "\n" );
afile.close();

outname='DATA__FREQ__AXIS__' + iparam + '.dat'
afile   = open( outname, 'w' );
for j in range(NtNew1):
    temp   =  str( '%.4e'%om1[j] )
    afile.write( temp );
    afile.write("\n");
afile.close();
plt.show();
