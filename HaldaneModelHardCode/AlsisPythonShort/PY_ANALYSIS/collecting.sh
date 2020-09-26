

vec=(0.000
0.063
0.126
0.188
0.251
0.314
0.377
0.440
0.503
0.565
0.628
0.691
0.754
0.817
0.880 
0.942
1.005
1.068
1.131
1.194
1.257
1.319
1.382
1.445
1.508
1.570
)

cd ./"SetData0"




cp ../AnalysisDichroism.py . 
chmod u+x AnalysisDichroism.py 


down=$1;
Nvec=$2;

####################################
###
for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

     echo ""
     echo "phi0 = ${vec[iparam]}"  
    
      "./"AnalysisDichroism.py "${vec[iparam]}" 220 

     echo ""
     echo "phi0 = ${vec[iparam]}"    

done


echo "done"


