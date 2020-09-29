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

bfile0 = "Dichroism__Phi0__"
bfile1 = "__M0t2__2.540__.dat"
bfile2 = "Haldane__Phase__Params__"

#Nlist  = 18
dichro = []
param  = []
index  = np.array([0,
                   0.06,
                   0.1,
                   0.2,
                   0.3,
                   0.4,
                   0.5,
                   0.6,
                   0.7,
                   0.8,
                   0.9,
                   1.00,
                   1.10,
                   1.16,
                   1.20,
                   1.30,
                   1.40,
                   1.50,
                   1.60,
                   1.70,
                   1.80,
                   1.90,
                   2.00,
                   2.10,
                   2.20,
                   2.30,
                   2.40,
                   2.50,
                   2.6,
                   2.7,
                   2.8,
                   2.9,
                   3.0,
                   3.1,
                   3.14
                   ])

print(len(index))
Nlist   = len(index)
M0t2    = 2.54

for n in range(Nlist):
    fname1 = bfile0 + str('%.3f'%index[n]) + bfile1
    fname2 = bfile2 + str('%.3f'%index[n]) + bfile1
    dichro.append( np.loadtxt(fname1) )
    param.append(  np.loadtxt(fname2) )


Norders = len(dichro[0][:,0])
print len(dichro[0][:,0])

cd          = np.zeros( (Nlist,Norders) )
li_ratio    = np.zeros( (Nlist,Norders) )
ri_ratio    = np.zeros( (Nlist,Norders) )
phase       = []
egap        = []
topoNumber   =[]
e0          = 0.0045;

for n in range(Nlist):
    phase.append(param[n][0])
    egap.append(param[n][3])
    topoNumber.append(param[n][2])
    for m in range(Norders):
        cd[n,m]=dichro[n][m,1]
        ri_ratio[n,m]=dichro[n][m,2]
        li_ratio[n,m]=dichro[n][m,3]


print  "phase", phase, "\n"
print "eg = ", egap,"\n"

######################################
#Ploting the current oscillations
width = 11
hight = width/1.62

fig     = plt.figure(figsize=(width,hight) )
ax1     = fig.add_axes([0.2, 0.15, 0.75, 0.75])
Nt      = len( cd[:,0] )

cd_av   = np.zeros((Nt,1))
Nv      = 8;
for i in range(3,Nv):
    cd_av[:,0] = cd_av[:,0]+cd[:,2*i]


#p1,=plt.plot(phase, cd[:,0], 'bo-', lw=2 );
#p2,=plt.plot(phase, cd[:,2], 'go-', lw=2 );
#p3,=plt.plot(phase, cd[:,4], 'ko-', lw=2 );
p4,=plt.plot(phase[1:Nlist], cd[1:Nlist,6],  'ko-', lw=2 );
p5,=plt.plot(phase[1:Nlist], cd[1:Nlist,8],  'yo-', lw=2 );
p6,=plt.plot(phase[1:Nlist], cd[1:Nlist,10], 'bo-', lw=4 );
p7,=plt.plot(phase[1:Nlist], cd[1:Nlist,12], 'go-', lw=4 );
p8,=plt.plot(phase[1:Nlist], cd[1:Nlist,14], 'mo-', lw=4 );
p9,=plt.plot(phase[1:Nlist], cd_av[1:Nlist]/(Nv-3), 'ro-', lw=8 );

plt.legend([p4, p5, p6, p7, p8], ['HO10','HO13','HO16','HO19','HO22'],loc='lower center');
plt.grid(True)

plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' Dichroism ', fontsize=26);
plt.tick_params(labelsize=28);

xmin =  0;
xmax =  np.pi;
ymin = -1.5;
ymax =  1.5;

plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );


ax1 = plt.gca()
ax2 = ax1.twiny()
pi=np.pi
    #for n in range( int(len(egap)/3) ):
#print egap[n*3]

ax2.set_xticks( np.arange(0,pi,pi/11) )
ax2.set_xticklabels(['','0.08','0.002','0.08','0.11','0.12','0.12','0.11','0.08','0.002','0.08'])
#ax1.set_xlabel('redshift')
ax2.set_xlim(xmin, xmax )

xmin =  0;
xmax =  np.pi;
ymin = -1.6;
ymax =  1.6;

plt.tick_params(labelsize=16);
plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );


fname='/Dichroism__2n+1__Mt2__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;
plt.savefig(fileName, dpi = 300);


######################################
######################################
fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

#p1,=plt.plot( phase, cd[:,1],  'bo-', lw=2 );
#p2,=plt.plot( phase, cd[:,3],  'go-', lw=2 );
#p3,=plt.plot( phase, cd[:,5],  'ko-', lw=2 );
#p4,=plt.plot( phase, cd[:,7],  'ro-', lw=2 );
#p5,=plt.plot( phase, cd[:,9],  'mo-', lw=2 );
p6,=plt.plot( phase, cd[:,11], 'yo-', lw=4 );
p7,=plt.plot( phase, cd[:,13], 'go-', lw=4 );


#plt.legend([p1, p2, p3,p4,p5,p6,p7], ['HO2', 'HO5', 'HO8','HO11','H14','HO17','HO20'],loc='lower right');

plt.legend([p6,p7], ['HO17','HO20'],loc='lower center');
plt.grid( True )



plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' Dichroism ', fontsize=26);

plt.tick_params(labelsize=28);
#plt.ylim( -1.05, 1.05 );

plt.ylim( ymin, ymax );
plt.xlim( xmin, xmax );


ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xticks( np.arange(0,pi,pi/11) )
ax2.set_xticklabels(['','0.08','0.002','0.08','0.11','0.12','0.12','0.11','0.08','0.002','0.08'])
ax2.set_xlim(xmin, xmax )




fname='/Dichroism__2n+2__Mt2__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;


plt.savefig(fileName, dpi = 300);


#################################
## Analysis Intensity ratio
width = 11
hight = width/1.62

fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

#p1,=plt.plot(phase, li_ratio[:,0],  'bo-', lw=2 );
#p2,=plt.plot(phase, li_ratio[:,2],  'go-', lw=2 );
p3,=plt.plot(phase, li_ratio[:,4],  'ko-', lw=2 );
p4,=plt.plot(phase, li_ratio[:,6],  'ro-', lw=2 );
p5,=plt.plot(phase, li_ratio[:,8],  'mo-', lw=2 );
p6,=plt.plot(phase, li_ratio[:,10], 'bo-', lw=4 );
p7,=plt.plot(phase, li_ratio[:,12], 'go-', lw=4 );
#p8,=plt.plot(phase, li_ratio[:,14], 'r.-', lw=4 );

#plt.legend([p1, p2, p3,p4,p5,p6,p7], ['HO1', 'HO4', 'HO7','HO10','HO13','HO16','HO19'],loc='lower right');

plt.legend([p3,p4,p5,p6,p7], [ 'HO7','HO10','HO13','HO16','HO19','HO22'],loc='lower center');
plt.grid(True)



plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' enhancement -left', fontsize=26);

plt.tick_params(labelsize=28);
#plt.ylim( 1e-08,1.e4 );
ax1.set_yscale('log')
plt.xlim( xmin, xmax )

ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xticks( np.arange(0,pi,pi/11) )
ax2.set_xticklabels(['','0.08','0.002','0.08','0.11','0.12','0.12','0.11','0.08','0.002','0.08'])
ax2.set_xlim(xmin, xmax )
plt.xlim( xmin, xmax );


fname='/Ratios_2n+1__Mt2__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;


plt.savefig(fileName, dpi = 300);


############################################
fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

#p1,=plt.plot(phase, ri_ratio[:,0],  'bo-', lw=2 );
#p2,=plt.plot(phase, ri_ratio[:,2],  'go-', lw=2 );
#p3,=plt.plot(phase, ri_ratio[:,4],  'ko-', lw=2 );
p4,=plt.plot(phase, ri_ratio[:,5],  'ro-', lw=2 );
p5,=plt.plot(phase, ri_ratio[:,7],  'mo-', lw=2 );
p6,=plt.plot(phase, ri_ratio[:,9], 'bo-', lw=4 );
p7,=plt.plot(phase, ri_ratio[:,11], 'go-', lw=4 );
#p8,=plt.plot(phase, ri_ratio[:,14], 'r.-', lw=4 );

plt.legend([p4,p5,p6,p7], ['HO08','HO11','HO14','HO17'],loc='lower center');
plt.grid(True);

plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' enhancement -right', fontsize=26);

plt.tick_params(labelsize=28);
#plt.ylim(1e-08,1.e4 );
ax1.set_yscale('log')
plt.xlim( xmin, xmax );
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xticks( np.arange(0,pi,pi/11) )
ax2.set_xticklabels(['','0.08','0.002','0.08','0.11','0.12','0.12','0.11','0.08','0.002','0.08'])
ax2.set_xlim(xmin, xmax )
plt.xlim( xmin, xmax );


fname='/Ratios_2n+2__Mt2__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;

plt.savefig(fileName, dpi = 300);


############################################
fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

p4,=plt.plot(phase, ri_ratio[:,5],  'ro-', lw=2 );
p5,=plt.plot(phase, ri_ratio[:,7],  'mo-', lw=2 );
p6,=plt.plot(phase, ri_ratio[:,9], 'bo-', lw=4 );
p7,=plt.plot(phase, ri_ratio[:,11], 'go-', lw=4 );


plt.legend([p4,p5,p6,p7], ['HO08','HO11','HO14','HO17'],loc='lower center');
plt.grid(True);

plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' enhancement -right', fontsize=26);

plt.tick_params(labelsize=28);
#plt.ylim(1e-08,1.e4 );
ax1.set_yscale('log')
plt.xlim( xmin, xmax );
ax1 = plt.gca()
ax2 = ax1.twiny()
ax2.set_xticks( np.arange(0,pi,pi/11) )
ax2.set_xticklabels(['','0.08','0.002','0.08','0.11','0.12','0.12','0.11','0.08','0.002','0.08'])
ax2.set_xlim(xmin, xmax )
plt.xlim( xmin, xmax );


fname='/Ratios_2n+2__Mt2__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;


plt.savefig(fileName, dpi = 300);



############################################
fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

p4,=plt.plot(phase, egap,  'ro-', lw=4 );
plt.grid(True);

plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' Eg (a.u.)', fontsize=26);

plt.tick_params(labelsize=28);
plt.xlim( xmin, xmax );

fname='/Energy-Gap__'+ str("%.3f"%M0t2)  + '.pdf'
fileName     = BasicPath + fname;
plt.savefig(fileName, dpi = 300);



############################################
fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])
p4,=plt.plot(phase, topoNumber,  'ro-', lw=4 );
plt.grid(True);

plt.xlabel('$\phi_0(rad)$ ', fontsize=26);
plt.ylabel(' Chern No.', fontsize=26);

plt.tick_params(labelsize=28);
plt.xlim( xmin, xmax );
plt.ylim(-0.2,1.2);

fname='/Chern__No__'+ str("%.3f"%M0t2)  + '.pdf';
fileName     = BasicPath + fname;

plt.savefig(fileName, dpi = 300);


plt.show();
