#!/bin/bash
down=$1
Nvec=$2


for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

     echo "${iparam}"
     
     
     "./"set_of_mpi_reduction.sh ${iparam} ${iparam} &

     echo ""
     

done
