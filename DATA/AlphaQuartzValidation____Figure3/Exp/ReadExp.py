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

plt.rc('axes',linewidth=3)

#Functions 
def masking_dipole(t,f,ta,tb,asigma,bsigma):
  Nt=len(t);
  amask=f;
  for j in range(0,Nt):
    if (t[j] <= ta):
      amask[j] = f[j]*np.exp(- ((t[j]-ta)/asigma)**2 );
      
    if (t[j] >= tb):
      amask[j] = f[j]*np.exp(- ((t[j]-tb)/bsigma)**2 );
      
  return amask;
  
#Getting parameters
set_of_params   = sys.argv;
lparam          = set_of_params[1]
#Creating path and file name for reading files 

BasicPath	    = os.getcwd();
FolderName      = "/" + lparam;
ProjPath  	    = BasicPath + FolderName;

xFileNameIntraInter  = '/' + 'Fig01.txt'

xInterIntraPath = BasicPath + xFileNameIntraInter

# LOADING DATA # 
xInterC         = np.loadtxt( xInterIntraPath );

Ntx             = len(xInterC[:,0])
#Nty = int(yInterC[0])
w0              = 1.5 ##Optical photon-freq. energy IR-LAser
horder          = xInterC[:,0]/w0

poli_total_inter1 = xInterC[:,1]
poli_total_inter2 = xInterC[:,2]
alphaq_total_inter = xInterC[:,3]


print( '\nNtx = ', Ntx)
#print '\nNty = ', Nty, '\n'

width = 11
hight = width/1.62
hzero = .971e-10

fig = plt.figure(figsize=(width,hight) )

ap = [0.163, 0.09, 0.823, 0.865]

ax1      = fig.add_axes(ap)#([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(3.2);



p1, = ax1.plot( horder, poli_total_inter1+hzero, color='b', lw=4, label='policrystal' );
plt.fill_between( horder, poli_total_inter1+hzero, 1e-0, color="skyblue", alpha=0.25 )


#p2, = plt.plot( horder, np.log10(alphaq_total_inter+hzero), 'y', lw=2, label='$\alpha$-quartz' );

plt.yscale('log')

xticks0  = np.arange(3,40,4);
plt.xticks( xticks0 );
plt.grid(True)

#plt.legend([p1, p2], [r'$\alpha$-quartz', 'policrystal'], fontsize=18);
plt.legend([p1], [r'$\alpha$-'+'quartz' ], fontsize=37);

fontz   = 37
xp0     = 3.5
yp0     = 0.4e5
ax1.text(xp0,yp0,'(a) Experiment',fontsize=fontz,color='k');

xaxmin      = 3;    # controling harmonic-order axis limits, down
xaxmax      = 25;   # controling harmonic-order axis limits, u
plt.xlim(xaxmin, xaxmax);
plt.ylim(0, 1e5);

plt.text(-.199, 0.5, r'$\rm I_{HHG}\,\,(arb.\,u.)$ ',
         fontsize=40,
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes);


plt.tick_params(labelsize=35);

fname = '/Figure3a.pdf'#'/HHG-Experimental-Spectrum' + '.pdf'
filename1           = fname;
fileNamePicture     = ProjPath  + filename1;
plt.savefig( fileNamePicture, dpi = 300 );


plt.show();

