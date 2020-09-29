#!/usr/bin/python
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib as mpl
import numpy as np
import os 
import sys 
import shutil 
import cmath
import errno
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d

plt.rc('axes',linewidth=3);
BasicPath            = os.getcwd();

bfile0 = "Dichroism__Phi0__1.571__M0t2__"
bfile1 = "__.dat"
bfile2 = "Haldane__Phase__Params__1.571__M0t2__"
dichro = []
param  = []
temp    = np.loadtxt("HPD__Param.dat");


M0t2    = temp[:,1];
index   = temp[:,1];
print( len(index) )
Nlist   = len(index);
for n in range(Nlist):
    fname1 = bfile0 + str('%.3f'%index[n]) + bfile1
    fname2 = bfile2 + str('%.3f'%index[n]) + bfile1
    dichro.append( np.loadtxt(fname1) )
    param.append(  np.loadtxt(fname2) )


Norders = len(dichro[0][:,0])
print( len(dichro[0][:,0]) )
cd          = np.zeros( (Nlist,Norders) )
li_ratio    = np.zeros( (Nlist,Norders) )
ri_ratio    = np.zeros( (Nlist,Norders) )
Mt2         = []
egap        = []
topoNumber  = []
e0          = 0.0045;

for n in range(Nlist):
    Mt2.append(param[n][1])
    egap.append(param[n][3])
    topoNumber.append(param[n][2])
    for m in range(Norders):
        cd[n,m]=dichro[n][m,1]
        ri_ratio[n,m]=dichro[n][m,2]
        li_ratio[n,m]=dichro[n][m,3]


print(  "Mt2", Mt2, "\n")
print( "eg = ", egap,"\n")

######################################
#Ploting the current oscillations
width = 12.
hight = width/1.62

fig     = plt.figure(figsize=(width,hight) )
#ax1     = fig.add_axes([0.075, 0.152, 0.89, 0.79])
ax      = fig.add_axes([0.1, 0.16, 0.89, 0.77]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.5);


Nt      = len( cd[:,0] )

cd_av   = np.zeros((Nt,1))
Nv      = 8;
for i in range(3,Nv):
    cd_av[:,0] = cd_av[:,0]+cd[:,2*i];

sw=8
int0=0


p4,=plt.plot(Mt2[int0:Nlist], cd[int0:Nlist,6],  'ks-', lw=2,markersize=sw );
p5,=plt.plot(Mt2[int0:Nlist], cd[int0:Nlist,8],  'ys-', lw=2,markersize=sw  );
p6,=plt.plot(Mt2[int0:Nlist], cd[int0:Nlist,10], 'bs-', lw=4,markersize=sw  );
p7,=plt.plot(Mt2[int0:Nlist], cd[int0:Nlist,12], 'gs-', lw=4,markersize=sw  );
#p8,=plt.plot(Mt2[int0:Nlist], cd[int0:Nlist,14], 'ms-', lw=4,markersize=sw  );
p9,=plt.plot(Mt2[int0:Nlist], cd_av[int0:Nlist]/(Nv-3), 'ro-', lw=5 );


plt.legend([p4, p5, p6, p7,  p9]
           , ['HO10','HO13','HO16','HO19','Average']
           ,loc=(.01,0.1),fontsize=22);


plt.xlabel(r'$ M_0/t_2$ ', fontsize=37);
plt.tick_params(labelsize=37);


xmin    = 0.;
xmax    = max(Mt2);
ymin    = -1.22;
ymax    = 1.45;

shift0  = 0.1;
xdown   = 0;
xup     = 3.*np.sqrt(3.);

plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );


plt.plot([xup,xup],     [ymin,ymax],    'r--',  lw=4);
plt.plot([xdown,xdown], [ymin,ymax],    'r--',  lw=4);


#plt.text(.2, 1.26, "(d)", fontsize=34, color="k");
plt.text(.25, 1.21, r'$\rm (d)\,\, CD\,\, vs \,\,M_0/t_2$',fontsize=36, color="k")


plt.text(1.3, 0.5,  r"$C=+1$", fontsize=27, color="k");
plt.text(7.5, 0.50, r"$C=0$",  fontsize=27, color="k");
ax1     = plt.gca();
ax2     = ax1.twiny();
pi      = np.pi;

Mt2New  = np.arange(0,max(Mt2),max(Mt2)/11.);
ax2.set_xticks( Mt2New );
Npt2    = len(Mt2New);

f       = interp1d( Mt2, egap, kind='cubic' );
egap1   = f(Mt2New);
ticks1  = [];

print( '\nLength egap = ', len(egap1), '   Np  = ', Npt2)


for n in range(0,Npt2):
    if n==0:
        ticks1.append(' ');
    else:
        ticks1.append('%.2f'%(egap1[n]*1.));


print( '\negap = ',egap1)

plt.yticks([-1,0,1]);
ax2.set_xticklabels(ticks1);
ax2.set_xlim( xmin, xmax );

shift0  = +1.0;
hd      = -1.1;

plt.fill_between( [xdown,xup-shift0], [1,1], [hd,hd], color=[0.0,0.8,0.99], alpha=0.04 );
plt.fill_between( [xup-shift0, xup+shift0], [1,1], [hd,hd], color=[0.8,0.0,0.], alpha=0.057 );
plt.fill_between( [xup+shift0, max(Mt2)], [1,1], [hd,hd], color=[0.0,0.9,0.], alpha=0.057 );

xmin    =  0.;
xmax    =  max(Mt2);
ymin    = -1.22;
ymax    =  1.45;#1.25;
plt.tick_params( labelsize=21 );
plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );
#fname='/Dichroism__2n+1__Phi__1.57__.pdf'

fname='/Figure5d.pdf'
fileName     = BasicPath + fname;
plt.savefig(fileName, dpi = 300);



plt.show();
