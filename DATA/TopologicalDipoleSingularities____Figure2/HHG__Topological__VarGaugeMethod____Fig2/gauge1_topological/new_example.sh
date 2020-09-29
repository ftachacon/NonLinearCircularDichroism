#!/bin/bash


####################################################
# This script runs antelope code for a basic example.
# You need to give a single "parameter" (directory name)
# via terminal and this scrip will compute:
# 1) evolution of the Semiconductor Bloch Eqs.
#     for Haldane Model
# 2) intra and inter band currents as a function of time



dirname=$1 # Pass a simple name to this scrip




##########################################
## You can modify the No. of MPI tasks

N_MPI_TASKs=23 #2   #No. of MPI task or cores for parallelization





###################################
mkdir ${dirname} #creating test directory


rm exec_hhg_mpi*


####################################
make #compiling code and generating executable 



####################################
#Copying exe* and other files to Test dir.
rm ./${dirname}/exec_hhg_mpi



cp exec_hhg_mpi ./${dirname}/
cp AlsisPythonShort/An* ./${dirname}/
cp AlsisPythonShort/drawSnapshots.py ./${dirname}/

cp input.dat ./${dirname}/


###################################
#Running test example
cd ./${dirname}


date
ls -lh exe*
time -p mpirun -n ${N_MPI_TASKs} exec_hhg_mpi

date 


#############################################
## Next line performe remain momentum integral
## along ky-direction #
##

#echo "--------------------"
#echo ""
#echo "Integration along ky direction via python script"

#time -p python postprocessingBForVis1.py / 0 ${N_MPI_TASKs}

#echo ""
#echo "End of the integration along ky-direction"


##removing temporary files
#rm full_*


###################################
##Analysis-first-step 
#cd  ./AlsisPythonShort/



######################################
##Convering data to python analysis codee
#./script_to_analysis.sh



cd ../
echo ""
echo "" 
