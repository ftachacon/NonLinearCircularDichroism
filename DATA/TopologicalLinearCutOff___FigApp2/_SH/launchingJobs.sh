#!/bin/bash


ellip=0.00       #Ellipticity--Laser--one
phi=$1           #Magnetic Haldane Model flux  


E0=$2            #0.0045       #Laser electric-field amplitude
T2=220           #Dephasing time 


name="EP__${ellip}__Phi0__${phi}_Efield__${E0}"   
mkdir ${name}
./copy_files.sh ./${name}


cd ./${name}


Mt2=2.54 #2.54
Nx=511 #401 #501  #691
Ny=595 #465 #291  #401
Nc=16   
dt=0.15


eps0=0.1
freg0=0
gauge=2

gbox_ky_down=0.01          #Botton Boundary  (lower) ky-BZ point for gauge modification
gbox_ky_up=0.70            #Upper  Boundary  (higher) ky-BZ point for gauge modification


fbz_shift=1                #Controlling BZ shift. it has to be 1 for yes, 0 for none 
ky_shift_down=-0.35        #Lowest BZ ky shift in a.u. 
factor_box=1

diagnostic=0

sbatch running-sbes2D.sh ${phi} ${Mt2} ${Nx} ${Ny} ${E0} ${Nc} ${dt} ${T2} ${ellip} ${eps0} ${gauge} ${freg0}  ${gbox_ky_down} ${gbox_ky_up} ${fbz_shift} ${ky_shift_down} ${diagnostic} ${factor_box} > id 

cd ../
