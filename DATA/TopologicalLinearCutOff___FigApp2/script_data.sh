#!/bin/bash
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
#################################
down=$1;
Nvec=$2;
for(( iparam=$down; iparam<=${Nvec}; iparam++ ))
do

    echo ""
    echo "dir = ${vec[iparam]}"

#    "./"AnalysisReader_PyDir.py 5 1 ${vec[iparam]} ${vec[iparam]} &
    "./"AnalysisReader_PyDir_Linear_Driver.py 5 1 ${vec[iparam]} ${vec[iparam]} &
    
    echo ""
    echo "dir = ${vec[iparam]}"

done
