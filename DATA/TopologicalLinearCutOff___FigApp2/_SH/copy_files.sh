#!/bin/bash 
module load python


dir=$1
#cp ./AlsisPythonShort/postprocessingBForVis1.py ./${dir} 
cp postprocessingBForVis1.py ./${dir}
cp running-sbes2D.sh ./${dir}
cp launch* ./${dir}
cp exec_hhg_mpi ./${dir}
cp copy_files.sh ./${dir}
cp *.sh ./${dir}


cd ./${dir}
ls -lh 

#python postprocessingBForVis1.py / 0 400 

#rm full_integrated_currents_rank_0* 

cd ..


echo "${dir}"
echo 
