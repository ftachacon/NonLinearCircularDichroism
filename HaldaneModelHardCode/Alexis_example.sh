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




#########################################
### SOME INPUTs PARAMETERS
## Set of Haldane M. Parameters

phi0=1.16           # Magnetic flux or phase of the complex 2nd hopping (rad.)


Mt2=2.54       # Ratio of on-site potential and 2nd hopping
gauge=2        # gauge parameter, can be 0, 1 and 2


rflag=0        # flag =0, No-reg.; flag=1, Taylor-Reg.; flag=2, Local-Gauss-Reg.
eps=0.1        # Reg. parameter for a fixed gauge, i.e. gauge=1

gbox_ky_down=0.01          #Botton Boundary  (lower) ky-BZ point for gauge modification
gbox_ky_up=0.70            #Upper  Boundary  (higher) ky-BZ point for gauge modification
ky_shift_down=-0.35        #Lowest BZ ky shift in a.u.


fbz_shift=1                 #Controlling BZ shift. it has to be 1 for yes, 0 for none
diag=1                      #diagnostic


Nx=250 #210 #401   #201 #201 #40#             # No. of points along kx
Ny=260 #240 #465   #233 #233 #117#            # No. of points along ky, ratio Nx/Ny=1.72 or, dpending on box 1.16



###############################
## Some laser field inputs ##
E0=0.0045 #0.0019           # Electric field strength (a.u.)
Ncycles=8         # No. of Opt. cycles
ellip=-1.0         # Ellipticity

dt=0.1             # Time-Steps (a.u.)
dephasing=220.      #1440. #220      # Dephasing (a.u.)

ksfactor=1
shotNumber=200


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



###################################
#Running test example
cd ./${dirname}


date
ls -lh exe*
time -p mpirun -n ${N_MPI_TASKs} exec_hhg_mpi ${phi0}  ${Mt2}  ${Nx}  ${Ny}  ${E0}  ${Ncycles}  ${dt}  ${dephasing} ${ellip} ${eps} ${gauge} ${rflag} ${gbox_ky_down} ${gbox_ky_up} ${fbz_shift} ${ky_shift_down} ${diag} ${ksfactor} ${shotNumber}

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
