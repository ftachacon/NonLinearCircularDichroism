#!/bin/bash

vec=(0.0
0.06
0.10
0.20
0.30
0.40
0.50
0.60
0.70
0.80
0.90
1.00
1.10
1.16
1.20
1.30
1.40
1.50
1.60
1.70
1.80
1.9
2.0
2.1
2.2
2.3
2.4
2.5
2.6
2.7
2.8
2.9
3.0
3.1
3.14
)

int=$1
Nf=$2


for(( iparam=${int}; iparam<=${Nf}; iparam++ ))
do

    echo ""
    echo "dir = ${vec[iparam]}"

    "./"AnalysisDichroism.py ${vec[iparam]} 220 

    echo ""
    echo "dir = ${vec[iparam]}"

done
