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

plt.rc('axes',linewidth=3);
BasicPath            = os.getcwd();


bfile0  = "Dichroism__Phi0__"
bfile1  = "__M0t2__2.540__T2__"
bfile2  = "Haldane__Phase__Params__"

#Nlist  = 18
dichro  = []
param   = []
phase0   = 0.06
T2   = np.array([220, 330, 440])

print(len(T2))
Nlist   = len(T2)
M0t2    = 2.54
for n in range(Nlist):
    fname1 = bfile0 + str('%.3f'%phase0) + bfile1 + str('%.3f'%T2[n]) + '__.dat'
    fname2 = bfile2 + str('%.3f'%phase0) + bfile1 + str('%.3f'%T2[n]) + '.dat'
    dichro.append( np.loadtxt(fname1) )
    param.append(  np.loadtxt(fname2) )


Norders = len(dichro[0][:,0])
print (len(dichro[0][:,0]) )
temp0=[]
cd_3n1  =[]
cd_3n2  =[]
order3n1=[]
order3n2=[]
phase       = []
egap        = []
topoNumber  = []
e0          = 0.0045;
No=int(Norders/2)
for i in range(Nlist):
    for n in range( No ):
        q=2*n
        q1=2*n+1
        cd_3n1.append( dichro[i][q,1] )
        order3n1.append( dichro[i][q,0] )
        order3n2.append( dichro[i][q,0] )
        cd_3n2.append( dichro[i][q1,1] )


print(  "phase", phase, "\n")
print( "eg = ", egap,"\n")
print( 'shape-orders = ', np.shape(order3n1) )

order3n1  = np.reshape(order3n1,(No,Nlist))
cd_3n1    = np.reshape(cd_3n1,(No,Nlist))

order3n2  = np.reshape(order3n2,(No,Nlist))
cd_3n2    = np.reshape(cd_3n2,(No,Nlist))
print( '\n', np.shape(order3n1) )

########################################
width = 11
hight = width/1.62
fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);


p1,=plt.plot(order3n1[:,0],cd_3n1[:,0],'o',color=[0,0,1], markersize=11,label=r'$T_2\,=\,220$')
p2,=plt.plot(order3n1[:,0],cd_3n1[:,1],'*',color=[0.0,1,0],     markersize=14,label=r'$\,\,\,= 330$',lw=8)
p3,=plt.plot(order3n1[:,0],cd_3n1[:,2],'s',color=[0.9,0,0],   markersize=9,label=r'$\,\,\,= 440$')

#plt.plot(order3n1,cd_3n1,'o',color=[1,0,0])

plt.legend([p1, p2, p3], [r'$T_2\,=\,\,220\,\rm a.u.$',r'$  \,\,\,\,=\,\, 330$',r'$  \,\,\,\,=\,\, 440$'],loc='upper right',fontsize=27);
#plt.grid(True)


#plt.xlabel(r'$\rm Harmonic--Order$ ', fontsize=34);
plt.ylabel(r'$\rm Circular\,\,\,Dichroism$', fontsize=34);



plt.grid()
xmin =  0;
xmax =  23;
ymin = -1.25;
ymax =  1.25;

xticks0  = np.arange(1,50,3);
plt.xticks( xticks0 );

plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );

plt.text(0.3, 1.05e-0, r'$\rm\,(c)\,Trivial\,\,3n+1$', fontsize=27, color="k");
#plt.text(1, 1e-0, r'$\rm\,\,\,\,\,\,\    Trivial\,\,\,3n+1$', fontsize=28, color="k");

plt.tick_params(labelsize=31);



fname='Figure7c.pdf'
#fileName     = BasicPath + fname;
plt.savefig(fname, dpi = 400);





############################################
width = 11
hight = width/1.62
fig     = plt.figure( figsize=(width,hight) );
ax      = fig.add_axes([0.16, 0.14, 0.81, 0.81]);
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3.0);


p1,=plt.plot(order3n2[:,0],cd_3n2[:,0],'o',color=[0,0,1],lw=2, markersize=11,label=r'$T_2\,=\,220$')
p2,=plt.plot(order3n2[:,0],cd_3n2[:,1],'*',color=[0.0,1,0],lw=14,     markersize=14,label=r'$  \,\,\,= 330$')
p3,=plt.plot(order3n2[:,0],cd_3n2[:,2],'s',color=[0.9,0,0],   markersize=9,label=r'$ \,\,\,= 440$')
plt.legend([p1, p2, p3], [r'$T_2\,=\,\,220\,\rm a.u.$',r'$  \,\,\,\,=\,\, 330$',r'$  \,\,\,\,=\,\, 440$'],loc='upper right',fontsize=27);



plt.grid();
xticks0  = np.arange(2,50,3);
plt.xticks( xticks0 );
plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );


plt.ylabel(r'$\rm Circular\,\,\,Dichroism$', fontsize=34);
plt.xlabel(r'$\rm Harmonic\,\,Order$', fontsize=34);
plt.text(0.3, 1.05e-0, r'$\rm (e)\,\,Trivial\,\,3n+2$', fontsize=27, color="k");
#plt.text(1, 1e-0, r'$\rm\,\,\,\,\,\,\    Trivial\,\,\,3n+2$', fontsize=28, color="k");


#fname='/Fig_CD_3n2_Trivial.pdf'
#fileName     = BasicPath + fname;

plt.tick_params(labelsize=31);

fname='Figure7e.pdf'
plt.savefig(fname, dpi = 400);





plt.show();
