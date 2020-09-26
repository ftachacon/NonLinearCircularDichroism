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
#import parameters as cs



##########################################
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
iparam 	    	=  float(set_of_params[1]);  #intensity param


#Creating path and file name for reading files 
BasicPath	        = os.getcwd();
FolderName          = "/" ;
ProjPath  	        = BasicPath + FolderName;

RPol        = 'HHG__Phi0__'+ str("%.3f"%iparam) + '__M0t2__2.540__e0__-1.000__.dat';
LPol        = 'HHG__Phi0__'+ str("%.3f"%iparam) +'__M0t2__2.540__e0__1.000__.dat';



#####################################################
FigureDir 	= '/' + 'Figures';
R_Path 	    = ProjPath  + RPol;
L_Path      = ProjPath  + LPol;


try:
  os.mkdir( ProjPath + FigureDir );
except OSError as exc:
  if (exc.errno != errno.EEXIST):
    raise exc; 
  pass  

try:
    os.mkdir( "SetData0" );
except OSError as exc:
    if (exc.errno != errno.EEXIST):
        raise exc;
    pass


#####################################################
# LOADING DATA # 
R_HHG 	        = np.loadtxt( R_Path );
L_HHG           = np.loadtxt( L_Path );


#####################################################
#shutil.copyfile( out, NewDir );
Nt          = int( L_HHG[0,0] );
Phi0        = L_HHG[0,1] ;
M0t2        = L_HHG[0,2] ;
ChernNo     = L_HHG[0,3] ;
Eg          = L_HHG[0,4] ;


#####################################################
print("\n+++++++++++++++++++++++++++++\nShort HALDANE M. STRUCTUTRE INFO.:" );
print("Eg        =   ", Eg, " a.u." );
print("Chern No. = ", ChernNo );
print("phi0      = ", Phi0, " rad." );
print("M0/t2     = ", M0t2 );
print('\n.......................++++++++++++++++++++++++=+');



#####################################################
### Frequency axis
wmax    = max(R_HHG[1:Nt,0]);
w     	= R_HHG[1:Nt,0];
dw      = w[1]-w[0];


#############################################################
## Computing harmonic spectra, inter and intra contribution
i 	    = cmath.sqrt(-1.);
hzero   = .98e-16;

R_DJ_m  = abs( R_HHG[1:Nt,1] + i * R_HHG[1:Nt,2] )**2 + hzero;
R_DJ_p  = abs( R_HHG[1:Nt,1] - i * R_HHG[1:Nt,2] )**2 + hzero;

L_DJ_m  = abs( L_HHG[1:Nt,1] + i * L_HHG[1:Nt,2] )**2 + hzero;
L_DJ_p  = abs( L_HHG[1:Nt,1] - i * L_HHG[1:Nt,2] )**2 + hzero;


############################
############################
##Up to 27th order
Nf      = 20;
CD      = np.zeros((Nf,1));
horder  = np.zeros((Nf,1));

horder[0] = 1
horder[1] = 2;
for n in range(2,Nf-1 ):
    horder[n]   = 3 + horder[n-2]
    horder[n+1] = 3 + horder[n-1]

print (horder)
#exit(0)

for n in range(0,Nf/2):
    
    #horder.append(3.*n)
    #horder[]=3*n+1;
    
    hint    =  horder[2*n,0] - 1./2.
    hfinal  =  horder[2*n+1,0] + 1./2.
    
    ind0    = np.floor( ( hint - w[0] )/dw );
    indf    = np.floor( ( hfinal - w[0] )/dw );
    #print ("i0 = ", ind0, "; if = ", indf)
    
    totalY     = np.sum( R_DJ_p[ind0:indf] + L_DJ_m[ind0:indf]  );
    asymmetryY = np.sum( R_DJ_p[ind0:indf] - L_DJ_m[ind0:indf]  )
    
    temp       = asymmetryY/totalY
    CD[n]    = temp;
#CD.append(temp)


#horder.shape = (Nf-1,1)
#CD.shape     = (Nf-1,1)

OutputData  = horder
OutputData  = np.concatenate( (OutputData, CD  ), axis=1 );


ofname          = "/Dichroism__Averange"  + "__Phi0__" + str("%.3f"%Phi0) + "__M0t2__" + str("%.3f"%M0t2) + "__.dat"

out                = ProjPath + ofname
np.savetxt(out, OutputData , fmt='%1.16e')

NewDir             = "./CollectingDC" + ofname
shutil.copyfile( out, NewDir );






#####################################################
############################
#Ploting the harmonic current radiations-oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])


p11, = plt.plot( w, np.log10( R_DJ_p ), 'b',  lw = 1.2, label='$D_tJ_{+,RCP}$' );
p22, = plt.plot( w, np.log10( L_DJ_m ), 'r',  lw = 1.2, label='$D_tJ_{-,RCP}$' );

plt.legend([p11, p22], [ '$D_tJ_{+,RCP}$', '$D_tJ_{-,RCP}$'],fontsize=18);


xp  = 1.3
yp  = -1.1


plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 );
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 );
plt.tick_params(labelsize = 28 );




xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

yticks0  = [-20., -15, -10.0, -5, 0, 5 ]#10., 15, 5.0];
plt.yticks( yticks0 );

xaxmin      = 0;    # controling harmonic-order axis limits, down
xaxmax      = 37;   # controling harmonic-order axis limits, u


plt.xlim( xaxmin, xaxmax );
plt.ylim( np.log10(hzero), 1.5 );
plt.grid(True)

fname='/LeftRightRotations' + "__Phi0__" + str("%.3f"%Phi0) + "__M0t2__" + str("%.3f"%M0t2)  + '.pdf'

filename     = FigureDir + fname;
fileName     = ProjPath  + filename;

plt.savefig(fileName, dpi = 300);




#####################################################
############################
#Ploting the harmonic current radiations-oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

plt.plot( horder, CD, 'go', lw = 2., label='$J_{x}$' );

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 );
plt.ylabel(r'$\rm CD $', fontsize=30 );
plt.tick_params(labelsize = 28 );


xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );

plt.grid(True)

xaxmin      = 0;    # controling harmonic-order axis limits, down
xaxmax      = 37;   # controling harmonic-order axis limits, u


plt.xlim(xaxmin, xaxmax);

fname='/CD_On_HHG_CNo' + "__Phi0__" + str("%.3f"%Phi0) + "__M0t2__" + str("%.3f"%M0t2)  + '.pdf'


filename     = FigureDir + fname;
fileName     = ProjPath  + filename;


plt.savefig(fileName, dpi = 300);


print( 'Eg = ',Eg, '\nphi0 = ', Phi0, '\nC = ',ChernNo)



plt.show();
#####################################################
