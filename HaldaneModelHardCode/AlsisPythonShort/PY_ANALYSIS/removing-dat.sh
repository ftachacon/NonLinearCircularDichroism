#!/bin/bash 

dir=$1
#cp postprocessingBForVis1.py ./${dir} 
#cp running-sbes2D.sh ./${dir}
#cp launch* ./${dir}
#cp exec_hhg_mpi ./${dir}



cd ./${dir}
ls -lh 

#python postprocessingBForVis1.py / 0 400 

#rm full_* 

#cp occupation__full__evol.dat ./AlsisPythonShort/

rm -R Figure 
rm connection.dat  #*.dat
rm curvature.dat 
rm dipoles.dat 
rm edispersion.dat 
rm gvelocities.dat 
rm mgrid.dat 
rm termial-outpput.dat

rm exec_hhg_mpi

cd ..


echo "${dir}"
echo 
