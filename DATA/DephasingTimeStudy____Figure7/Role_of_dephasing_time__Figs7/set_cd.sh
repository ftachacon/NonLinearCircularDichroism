#!/bin/bash

vec=(220
330
440
)

int=0
Nf=2

phi=$1

for(( iparam=${int}; iparam<=${Nf}; iparam++ ))
do

    echo ""
    echo "dir = ${vec[iparam]}"

    "./"AnalysisDichroism.py ${phi} ${vec[iparam]}

    echo ""
    echo "dir = ${vec[iparam]}"

done
