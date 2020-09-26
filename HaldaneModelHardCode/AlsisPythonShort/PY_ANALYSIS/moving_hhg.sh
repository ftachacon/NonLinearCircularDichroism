#!/bin/bash


T2=$1
fname="Dephasing_T2_${T2}"

mkdir ${fname}

f1="Efield0.005___Phi1.16___Dep${T2}___Ellip-01"
f2="Efield0.005___Phi1.16___Dep${T2}___Ellip+01"


cp ./${f1}/HHG* ./${fname}
cp ./${f2}/HHG* ./${fname}

