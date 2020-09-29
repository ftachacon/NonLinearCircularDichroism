
mpitasks=$2
f=$1 

module load python

cp postprocessing* ./$f
cd ./$f

ls -lh full* 

time -p ./postprocessingBForVis1.py / 0 ${mpitasks}

ls -lh int*

echo ""
echo " Reduction is finished "
echo " "

pwd 

cd ../
