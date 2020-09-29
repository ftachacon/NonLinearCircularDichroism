#!/bin/bash
ellip=$1
phase=$2
Mt2=$3

fdir="EP__${ellip}__Phi__${phase}__Mt2__${Mt2}"
mkdir ${fdir}

cp -f exec* ${fdir}
cp -f mpi_gcc-4.9.2.sh ${fdir}
cp -f input.dat ${fdir}


