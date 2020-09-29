#!/bin/bash
vec=(
Efield0.005___Phi0.06___Dep110___Ellip+01
Efield0.005___Phi0.06___Dep110___Ellip-01
Efield0.005___Phi0.06___Dep220___Ellip+01
Efield0.005___Phi0.06___Dep220___Ellip-01
Efield0.005___Phi0.06___Dep330___Ellip+01
Efield0.005___Phi0.06___Dep330___Ellip-01
Efield0.005___Phi0.06___Dep440___Ellip+01
Efield0.005___Phi0.06___Dep440___Ellip-01
Efield0.005___Phi0.06___Dep880___Ellip+01
Efield0.005___Phi0.06___Dep880___Ellip-01
Efield0.005___Phi1.16___Dep220___Ellip+01
Efield0.005___Phi1.16___Dep220___Ellip-01
Efield0.005___Phi1.16___Dep330___Ellip+01
Efield0.005___Phi1.16___Dep330___Ellip-01
Efield0.005___Phi1.16___Dep440___Ellip+01
Efield0.005___Phi1.16___Dep440___Ellip-01
Efield0.005___Phi1.16___Dep880___Ellip+01
Efield0.005___Phi1.16___Dep880___Ellip-01
)

down=$1
Nvec=$2

for(( iparam=${down}; iparam<=${Nvec}; iparam++ ))
do

    echo ""
    echo "dir = ${vec[iparam]}"

    "./"AnalysisReader_PyDir.py 5 01 ${vec[iparam]} ${vec[iparam]}

    echo "i = ${iparam}"
    echo "dir = ${vec[iparam]}"

done

