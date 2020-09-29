#!/bin/bash

           
Mr=2.54    #Local potential 
Nx=401     #No. of kx points 
Ny=151   
e0=0.006   #laser field strength
Ncy=10     #No. of Opt. cycles
dt=1.25    #time-step 

down=$1     # initial-phase 
upper=$2    # final-phase 


ex=exec_ac_hhg
phi0=3 #-3.141592654
dp=0.125663706
run=running-sbes2D.sh


#LOOP ON MAGNETIC FLUX \phi
for((iparam=${down}; iparam<=${upper}; iparam++ ))
do

    dir="Phi${iparam}"

    mkdir ${dir}
    cp ${ex} ${dir}/
    cp ${run} ${dir}
    cd ${dir}/
    
    ls -lh 

    echo ""
    echo ""
    echo "Phi index = ${iparam}"
    echo "Hello from dir. ${dir}"
    echo "Mr = $Mr;  Nx = ${Nx};  Ny = ${Ny};   e0 = ${e0};  Ncy = ${Ncy};  dt = ${dt}\n"
    echo ""
    
     sbatch running-sbes2D.sh ${iparam} ${Mr} ${Nx} ${Ny} ${e0} ${Ncy} ${dt} > id &
     
     echo "job id"
     more id 
     echo ""
     cd ../
     
    #cp down_load* $num_sim/

done
