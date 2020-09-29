#!/bin/bash

vec=(EP__+01__Phi__0.20__Mt2__2.54
EP__+01__Phi__0.30__Mt2__2.54
EP__+01__Phi__0.40__Mt2__2.54
EP__+01__Phi__0.50__Mt2__2.54
EP__+01__Phi__0.60__Mt2__2.54
EP__+01__Phi__0.70__Mt2__2.54
EP__+01__Phi__1.00__Mt2__2.54
EP__+01__Phi__2.70__Mt2__2.54
EP__-01__Phi__0.20__Mt2__2.54
EP__-01__Phi__0.30__Mt2__2.54
EP__-01__Phi__0.40__Mt2__2.54
EP__-01__Phi__0.50__Mt2__2.54
EP__-01__Phi__0.60__Mt2__2.54
EP__-01__Phi__0.70__Mt2__2.54
EP__-01__Phi__1.00__Mt2__2.54
EP__-01__Phi__2.70__Mt2__2.54
EP__00__Phi__0.00__Mt2__2.54
EP__00__Phi__1.57__Mt2__2.54
)


for(( iparam=$1; iparam<=$2; iparam++ ))
do

    echo ""
    echo "dir = ${vec[iparam]}"

    rm "./${vec[iparam]}/"antelope*
    rm "./${vec[iparam]}/"e*
    rm -r "./${vec[iparam]}/"NPhi*
    rm -r "./${vec[iparam]}/"P*
    rm -r "./${vec[iparam]}/"c*
    rm -r "./${vec[iparam]}/"d*
    rm -r "./${vec[iparam]}/"g*
    rm -r "./${vec[iparam]}/"HHG*
    rm -r "./${vec[iparam]}/"stdout.dat
    rm -r "./${vec[iparam]}/"mgrid.dat

    echo ""
#echo "dir = ${vec[iparam]}"

done

