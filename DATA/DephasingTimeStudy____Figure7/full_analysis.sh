
vec=(
EP__+01__Phi0__0.00
EP__+01__Phi0__0.06
EP__+01__Phi0__0.10
EP__+01__Phi0__0.20
EP__+01__Phi0__0.30
EP__+01__Phi0__0.40
EP__+01__Phi0__0.50
EP__+01__Phi0__0.60
EP__+01__Phi0__0.80
EP__+01__Phi0__0.90
EP__+01__Phi0__1.10
EP__+01__Phi0__1.16
EP__+01__Phi0__1.20
EP__+01__Phi0__1.30
EP__+01__Phi0__1.40
EP__+01__Phi0__1.50
EP__+01__Phi0__1.60
EP__+01__Phi0__1.70
EP__+01__Phi0__1.80
EP__+01__Phi0__1.90
EP__+01__Phi0__2.00
EP__+01__Phi0__2.10
EP__+01__Phi0__2.20
EP__+01__Phi0__2.30
EP__+01__Phi0__2.40
EP__+01__Phi0__2.50
EP__+01__Phi0__2.60
EP__+01__Phi0__2.80
EP__+01__Phi0__2.90
EP__+01__Phi0__3.00
EP__+01__Phi0__3.10
EP__+01__Phi0__3.14
EP__-01__Phi0__0.00
EP__-01__Phi0__0.06
EP__-01__Phi0__0.10
EP__-01__Phi0__0.20
EP__-01__Phi0__0.30
EP__-01__Phi0__0.40
EP__-01__Phi0__0.50
EP__-01__Phi0__0.60
EP__-01__Phi0__0.80
EP__-01__Phi0__0.90
EP__-01__Phi0__1.10
EP__-01__Phi0__1.16
EP__-01__Phi0__1.20
EP__-01__Phi0__1.30
EP__-01__Phi0__1.40
EP__-01__Phi0__1.50
EP__-01__Phi0__1.60
EP__-01__Phi0__1.70
EP__-01__Phi0__1.80
EP__-01__Phi0__1.90
EP__-01__Phi0__2.00
EP__-01__Phi0__2.10
EP__-01__Phi0__2.20
EP__-01__Phi0__2.30
EP__-01__Phi0__2.40
EP__-01__Phi0__2.50
EP__-01__Phi0__2.60
EP__-01__Phi0__2.80
EP__-01__Phi0__2.90
EP__-01__Phi0__3.00
EP__-01__Phi0__3.10
EP__-01__Phi0__3.14
)



####################################
###

down=$1
Nvec=$2 

for(( iparam=$down; iparam<=${Nvec}; iparam++ ))
do 

    echo ""
    echo "dir = ${vec[iparam]}"  
     
    "./"AnalysisReader_PyDir.py 5 01 ${vec[iparam]} ${vec[iparam]} 
     
    echo ""
    echo "dir = ${vec[iparam]}"    

done


