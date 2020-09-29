#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:53:15 2019

@author: achacon
"""

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
from scipy.integrate import odeint
from scipy.integrate import ode
import time
from scipy import special
from scipy.integrate import quadrature, quad
#import parameters as cs


################################
plt.rc('axes',lw=2.2)


## CONSTANTS
pi      = np.pi;



################################
### INPUT PARAMETS
set_of_params   = sys.argv;
iparam          = set_of_params[1]  #directory name
kparam          = int(set_of_params[2]); #down index
nparam          = int(set_of_params[3]); #upper-index



################################
##  Basic current path
BasicPath       = os.getcwd();



################################
#Need to load files
FolderName              = "/"+iparam;
FolderProjPath          = BasicPath + FolderName;

flname              = "outlaserdata.dat";
pathl               = FolderProjPath + flname;
laser               = np.loadtxt( pathl );
taxis 	 	    = laser[:,0];
Ntime 		    = len(taxis)

#####################################################
## Creating figure dir.
#FigureDir               = '/Figure';
#try:
#    os.mkdir(BasicPath+FigureDir)#0777);
#except OSError as exc:
#    if (exc.errno != errno.EEXIST):
#        raise exc;
#        pass


inter_dip_ky 	= np.zeros( (Ntime,2) )
intra_jc_ky  	= np.zeros( (Ntime,2) )
occup0       	= np.zeros( (Ntime,2) )
occup0[:,0]  	= laser[:,0]



set          = range(kparam,nparam)
print( "Set of ky-index = ", set, "\nNtime = ", Ntime );






########################################
for n in set:

    fname0  = 'full_integrated_currents_rank_' + str('%.6d' % n) + ".dat"  #+  '_kyIndex_' + str('%.6d' % n) + '.dat'
    path0       = FolderProjPath + fname0;    
         

    ########################################
    ### LOADING DATA #################
    inter_intra_cts =  np.loadtxt( path0 );
    

    
    inter_dip_ky[:,0]+= inter_intra_cts[:,1];
    inter_dip_ky[:,1]+= inter_intra_cts[:,2];

    intra_jc_ky[:,0]+=  inter_intra_cts[:,3];
    intra_jc_ky[:,1]+=  inter_intra_cts[:,4];

    occup0[:,1]+= inter_intra_cts[:,5]



oname0      = "/interband_dipole_full_evol.dat";
oname1      = "/intraband_current_full_evol.dat";
oname2      = "/occupation__full__evol.dat";



qpath0      = BasicPath + oname0;
qpath1      = BasicPath + oname1;
qpath2      = BasicPath + oname2;


np.savetxt( qpath0, inter_dip_ky, '%.16e' );
np.savetxt( qpath1, intra_jc_ky, '%.16e' );
np.savetxt( qpath2, occup0, '%.16e' );


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
########################################

