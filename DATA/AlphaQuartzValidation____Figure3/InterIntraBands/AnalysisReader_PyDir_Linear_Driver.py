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



#Clear
#plt.cla()
#plt.clf()

#Getting parameters
set_of_params   = sys.argv;

iparam 	    	=  int( set_of_params[1] );  #intensity param
kparam          =  int( set_of_params[2] );  # Chern-No.
lparam      	=  set_of_params[3] ;        #dir-name
mparam      	=  set_of_params[4] ;        #xdirection





#Creating path and file name for reading files 
BasicPath	        = os.getcwd();

FolderName              = "/" + lparam;
ProjPath  	        = BasicPath + FolderName;


FileNameIntraInter      = '/interband_dipole_full_evol.dat';
FileNameIntraInter1     = '/intraband_current_full_evol.dat';
FileNameLaser 	        = '/outlaserdata.dat';
FileParameters 		    = '/setOfparameters.dat';
set_DataName 	        = '/SetData0';



#####################################################
FigureDir 	            = '/' + mparam + 'Figure';
IntraInterPath 	        = ProjPath  + FileNameIntraInter;
IntraInterPath1         = ProjPath  + FileNameIntraInter1;
LaserPath 	            = ProjPath  + FileNameLaser;
ParamPath 	            = ProjPath  + FileParameters;



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


#########################
# LOADING DATA # 
InterC 	        = np.loadtxt( IntraInterPath );
IntraC          = np.loadtxt( IntraInterPath1 );
Laser 	        = np.loadtxt( LaserPath );
Params 		    = np.loadtxt( ParamPath );



#####################################################
#shutil.copyfile( out, NewDir );
Phi0        = Params[2,0] ;
M0          = Params[2,1] ;
t1          = Params[1,5] ;
t2          = Params[1,6] ;
M0t2        = M0/t2;
ChernNo     = Params[2,4] ;
Eg          = Params[2,2]; 
ellip       = Params[3,1];

if abs(ChernNo) < 1.e-4:
    ChernNo=0.;


#####################################################
print("\n+++++++++++++++++++++++++++++\nHALDANE M. STRUCTUTRE INFO.:")
print("Eg        = ", Eg, " a.u."   );
print("Chern No. = ", ChernNo       );
print("phi0      = ", Phi0, " rad." );
print("M0        = ", M0    );
print("t1        = ", t1    );
print("t2        = ", t2    );
print("M0/t2     = ", M0t2  );


print("\n++++++++++++++++++++++++++++\nLaser features:" );
print("E0     = ", Params[0,0], " a.u." );
print("I0     = ", Params[0,1], " W/cm^2" );
print("w0     = ", Params[0,2], " a.u." );
print("N.O.C. = ", Params[0,3], " " );
print("Ellip  = ", ellip, " " );
print("\n\n Time-step dt  = ", Params[0,6], " a.u. " );
print("Total No. of time-steps  = ", Params[0,5], " \n" );



print('\n+=++++++++++++++++++++++++=+');
#####################################################


#print "size dim of ByHand", ByHand.shape
#Arranging/organizing data and parameters 
#Nt			= cs.Ntime;



#frequency and laser period
w0		= Params[0,2];
T0		= 2.*np.pi/w0;



print( "period: ", T0, " a.u.")
print( "mean-freq: ", w0, " a.u.")


#####################################################
#Laser intensity
E0          = np.sqrt( Params[0,1]/3.5e16);


#####################################################
#Time axis and electric field 
t		    = Laser[:,0];
Efield      = Laser[:,1];

Nt          = len(t);
dt 	        = Params[0,6];#t[1]-t[0];


print( "dt: ", dt, " a.u." );
print( "Ntime: ", Nt );


#####################################################
#Frequency axis 
wmax 	    	= np.pi/dt;
dw		        = 2.*wmax/float(Nt);
wmin  	    	= -wmax;
wmax 	    	= +wmax;
w		= np.arange( wmin, wmax-dw, dw );
w           	= np.linspace( wmin, wmax-dw, num = Nt )




#####################################################
if (Nt != len(w)):
    print ("\nsNtime have to be equal to length (w)\n")
    exit()
#####################################################




print( "Lengths of omega axis, Nomega: ", w.shape)
print( "dw= ", dw , " a.u.")
print( '\nwmax={0:.2f}'.format(wmax) );


#####################################################



#####################################################
print( "\nlen of interC: ", InterC.shape)
print('\n---')




#####################################################
#Momentum Integration for the "remains" of MPI code
xinterC 	= InterC[:,0] - InterC[0,0];
yinterC     = InterC[:,1] - InterC[0,1];

xJintra  	= IntraC[:,0] - IntraC[0,0];
yJintra     = IntraC[:,1] - IntraC[0,1];




#####################################################
#Computing inter-current contribution from dipole polarization InterC1
xJinter		   = np.zeros( xinterC.shape, np.float );
yJinter        = np.zeros( yinterC.shape, np.float );



xJinter 	   = xinterC    #np.diff( xinterC )/dt;
yJinter        = yinterC    #np.diff( yinterC )/dt;


#
######################################################
##Ploting the current oscillations
#width = 11
#hight = width/1.62
#
#
#fig = plt.figure(figsize=(width,hight) )
#
#
#p1, = plt.plot( t, Efield/max(Efield)*max(xJinter), 'red', lw=2 );
#p2, = plt.plot( t, xJinter, 'green', lw = 1.5 );
#p3, = plt.plot( t, yJinter, 'blue',  lw = 1.5 );
#plt.legend([p1, p2, p3], ['$E_{L}$', '$J_{x,er}$', '$J_{y,er}$']);
#plt.xlabel('time (a.u.) ',  fontsize=18);
#plt.ylabel('Jx/Jy (a.u.) ', fontsize=18);
#plt.tick_params(labelsize=18);
#
#xaxmin      = t.min();    #
#xaxmax      = t.max();    #
#plt.xlim( xaxmin, xaxmax );
#plt.tight_layout();
#filename0           = FigureDir + '/CurrentOscillations.pdf';
#fileNamePicture     = ProjPath + filename0;
#plt.savefig( fileNamePicture );
#





#####################################################
#Populations
width = 11
hight = width/1.62


#############################################################
## Computing harmonic spectra, inter and intra contribution
i 	= cmath.sqrt(-1)#complex(0,1);
print ("\ncomplex number: ", i)


#############################################################
#filtering dipole and current oscillations by means of 
#applying a smoth time mask over the beginning and end of pulses, this will avoid 
#high esporeous frequencies...
ta 	    = t[0]  + T0*18.0;#T0*3.0;
tb 	    = t[-1] - T0*18.0;#T0*3.0;

asigma  = T0*4
bsigma  = T0*4


xJinterMasked = np.blackman(Nt)*xJinter #
yJinterMasked = np.blackman(Nt)*yJinter #


xJintraMasked = np.blackman(Nt)*xJintra #
yJintraMasked = np.blackman(Nt)*yJintra #



#####################################################
##Calculating FFT 
xFFT_Jinter 	= fft( xJinterMasked )*dt;
yFFT_Jinter 	= fft( yJinterMasked )*dt;

xFFT_Jintra     = fft( xJintraMasked )*dt;
yFFT_Jintra     = fft( yJintraMasked )*dt;

xFFT_Jinter      = -np.fft.fftshift( xFFT_Jinter )*(1j)*w;
yFFT_Jinter      = -np.fft.fftshift( yFFT_Jinter )*(1j)*w;

xFFT_Jintra  	 = -np.fft.fftshift( xFFT_Jintra )*(1j)*w;
yFFT_Jintra  	 = -np.fft.fftshift( yFFT_Jintra )*(1j)*w;

xFullRadiation 	 = xFFT_Jinter   +  xFFT_Jintra;
yFullRadiation   = yFFT_Jinter   +  yFFT_Jintra;

dJ_p             = xFullRadiation - (1j)*yFullRadiation
dJ_m             = xFullRadiation + (1j)*yFullRadiation


############################
############################
hzero           = .98e-16

#xSinter         = np.log10( abs( xFFT_Jinter )**2 + hzero );
#ySinter         = np.log10( abs( yFFT_Jinter )**2 + hzero );
#Sinter          = np.log10( abs( xFFT_Jinter )**2 + abs( yFFT_Jinter )**2 + hzero );


xSinter           = abs( xFFT_Jinter )**2 + hzero ;
ySinter           = abs( yFFT_Jinter )**2 + hzero ;
xTotalInterIntra = abs( xFullRadiation )**2 + hzero;
Sinter            = abs( xFFT_Jinter )**2 + abs( yFFT_Jinter )**2 + hzero ;


#xSintra         = np.log10( abs( xFFT_Jintra )**2 + hzero );
#ySintra         = np.log10( abs( yFFT_Jintra )**2 + hzero );
#Sintra          = np.log10( abs( xFFT_Jintra )**2 + abs( yFFT_Jintra )**2 + hzero);

xSintra             =  abs( xFFT_Jintra )**2 + hzero
ySintra             =  abs( yFFT_Jintra )**2 + hzero ;
yTotalInterIntra    =  abs( yFFT_Jinter )**2 + abs( yFFT_Jintra )**2 + hzero;
Sintra              =  abs( xFFT_Jintra )**2 + abs( yFFT_Jintra )**2 + hzero;



xSpectrum 	    = np.log10( abs( xFullRadiation )**2 + hzero );
ySpectrum       = np.log10( abs( yFullRadiation )**2 + hzero );
Spectrum        = np.log10( abs( xFullRadiation )**2 + abs( yFullRadiation )**2 + hzero );


a_dJ_p          = np.log10( abs( dJ_p )**2 + hzero );
a_dJ_m          = np.log10( abs( dJ_m )**2 + hzero );
#####################################################




print( "Spectra shape= ",  xSinter.shape )
print( ";     frequency axis shape= ", w.shape )


w.shape             = (Nt,1)
xSpectrum.shape     = (Nt,1)
ySpectrum.shape     = (Nt,1)
Spectrum.shape      = (Nt,1)

xFullRadiation.shape = (Nt,1);
yFullRadiation.shape = (Nt,1);



Nthalf          = int(np.floor(Nt/2));
temp0           = np.array([Nthalf]);
temp0.shape     = (1,1)
w_output        = np.concatenate( ( temp0, w[Nthalf:Nt-1]/w0 ),axis=0 );


temp0           = np.array([Phi0]);
temp0.shape     = (1,1)
real_jxoutput   = np.concatenate( ( temp0, np.real(xFullRadiation[Nthalf:Nt-1]) ),axis=0 );


temp0           = np.array([M0t2]);
temp0.shape     = (1,1)
imag_jxoutput    = np.concatenate( ( temp0, np.imag(xFullRadiation[Nthalf:Nt-1]) ),axis=0 );



temp0           = np.array([ChernNo]);
temp0.shape     = (1,1)
real_jyoutput   = np.concatenate( ( temp0, np.real( yFullRadiation[Nthalf:Nt-1]) ),axis=0 );


temp0           = np.array([Eg]);
temp0.shape     = (1,1)
imag_jyoutput   = np.concatenate( ( temp0, np.imag( yFullRadiation[Nthalf:Nt-1] ) ),axis=0 );


temp0           = np.array([dt]);
temp0.shape     = (1,1)
spectrum_output = np.concatenate( ( temp0, pow( 10., Spectrum[Nthalf:Nt-1] ) ),axis=0 );


OutputData      = w_output;
OutputData      = np.concatenate( (OutputData, real_jxoutput   ), axis=1 );
OutputData      = np.concatenate( (OutputData, imag_jxoutput   ), axis=1 );
OutputData      = np.concatenate( (OutputData, real_jyoutput   ), axis=1 );
OutputData      = np.concatenate( (OutputData, imag_jyoutput   ), axis=1 );
OutputData      = np.concatenate( (OutputData, spectrum_output ), axis=1 );


ofname          = "/HHG__Phi0__" + str("%.3f"%Phi0) + "__M0t2__" + str("%.3f"%M0t2) + "__e0__" + str("%.3f"%ellip) +"__E0" +str("%.4f"%E0) + "__.dat"

out                = ProjPath + ofname
np.savetxt(out, OutputData , fmt='%1.16e')


NewDir = "./SetData0" + ofname
shutil.copyfile( out, NewDir );


#####################################################
############################
##FIGURE 1
#Ploting the harmonic current radiations-oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
#ap = [0.1865, 0.1511, 0.811, 0.814]
ap = [0.1865, 0.152, 0.81, 0.81]
ax0     = plt.axes(ap);

for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(3);


p1, = plt.plot( w/w0, xSintra, 'r',  lw = 3 );
p2, = plt.plot( w/w0, xSinter, '-b',  lw = 2.5 );
p3, = plt.plot( w/w0, xTotalInterIntra, '--g',  lw = 3 );

plt.legend([p1, p2, p3], [ r'$\rm I_{ra}^{(x)}$', r'$\rm I_{er}^{(x)}$', r'$\rm I_{total}^{(x)}$'],fontsize=31);

plt.yscale('log');

plt.xlabel(r'$\rm Harmonic\,Order$', fontsize=34 );


plt.text(-.229, 0.5, r'$\rm I_{HHG}\,\,(arb.\,u.)$ ',
         fontsize=39,
         horizontalalignment='left',
         verticalalignment='center',
         rotation=90,
         clip_on=False,
         transform=plt.gca().transAxes)

plt.tick_params(labelsize = 35 );
xticks0  = np.arange(1,50,4);
plt.xticks( xticks0 );


yticks0  = [1e-16, 1e-12, 1e-8, 1e-4, 1e0 ];
plt.yticks( yticks0 );


xaxmin      = 0;    # controling harmonic-order axis limits, down
xaxmax      =26.5;   # controling harmonic-order axis limits, u
plt.xlim( xaxmin, xaxmax );
plt.ylim([0.98e-16,1.00e-0])
plt.grid(True)

plt.text(0.65, 3.50e-2, '(c) Par. Emission ', fontsize=34, color="k");


fname='/Figure3c.pdf'

filename1           = FigureDir + fname;
fileNamePicture     = ProjPath  + filename1; 

plt.savefig(fileNamePicture, dpi = 300);
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname );





#####################################################
############################
#Ploting the harmonic current radiations-oscillations
width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ap = [0.1865, 0.147, 0.811, 0.814]
ax1     = plt.axes(ap);

for axis1 in ['top','bottom','left','right']:
    ax1.spines[axis1].set_linewidth(3);


p1, = plt.plot( w/w0, ySintra, 'r',  lw = 3 );
p2, = plt.plot( w/w0, ySinter, '-b',  lw = 2.5 );
p3, = plt.plot( w/w0, yTotalInterIntra, '--g',  lw = 3 );

plt.legend([p1, p2, p3], [ r'$\rm I_{ra}^{(y)}$', r'$\rm I_{er}^{(y)}$', r'$\rm I_{total}^{(y)}$'],fontsize=31);


plt.yscale('log')
plt.xlabel(r'$\rm Harmonic\,Order$', fontsize=34 );
plt.tick_params(labelsize = 35 );

plt.yscale('log')


xaxmin      = 0;    # controling harmonic-order axis limits, down
xaxmax      =26.5;   # controling harmonic-order axis limits, us
plt.xlim( xaxmin, xaxmax );
xticks0  = np.arange(1,26,4);
plt.xticks( xticks0 );
plt.grid(True)


yticks0  = [1e-16, 1e-12, 1e-8, 1e-4, 1e0 ];
plt.yticks( yticks0 );
plt.ylim([0.98e-16,1.00e-0])
plt.text(0.65, 3.50e-2, '(d) Per. Emission ', fontsize=34, color="k");


fname='/Figure3d.pdf'
filename1           = FigureDir + fname;
fileNamePicture     = ProjPath  + filename1; 


plt.savefig(fileNamePicture, dpi = 300);
shutil.copyfile( fileNamePicture, BasicPath + set_DataName + fname );


print( 'Eg = ',Eg, '\nM0 = ',M0, '\nphi0 = ', Phi0, '\nC = ',ChernNo, '\nmparam =', mparam, '-direction')


#
plt.show();
#####################################################
