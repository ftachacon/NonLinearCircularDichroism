#!/bin/bash

module load python 

vec=(EP__0.00__Phi0__0.06_Efield__0.002
EP__0.00__Phi0__0.06_Efield__0.003
EP__0.00__Phi0__0.06_Efield__0.004
EP__0.00__Phi0__0.06_Efield__0.0045
EP__0.00__Phi0__0.06_Efield__0.005
EP__0.00__Phi0__0.06_Efield__0.006
EP__0.00__Phi0__0.06_Efield__0.007
EP__0.00__Phi0__1.16_Efield__0.002
EP__0.00__Phi0__1.16_Efield__0.003
EP__0.00__Phi0__1.16_Efield__0.004
EP__0.00__Phi0__1.16_Efield__0.0045
EP__0.00__Phi0__1.16_Efield__0.005
EP__0.00__Phi0__1.16_Efield__0.006
EP__0.00__Phi0__1.16_Efield__0.007
)


####################################
###

mpitasks=119 
Nvec=$2 
down=$1
file="intraband_current_full_evol.dat"


####################################
###
for((iparam=$down; iparam<=${Nvec}; iparam++ ))
do 

     echo ""
     echo "dir = ${vec[iparam]}"  
    
    "./"reduction.sh  "${vec[iparam]}"  ${mpitasks}
    "./"removing_mpi_files.sh  "${vec[iparam]}"   

     echo ""
     echo "dir = ${vec[iparam]}"    

done


echo ""
echo "\nEnd of the launching phase study at ${vec[0]}\n"
echo ""

