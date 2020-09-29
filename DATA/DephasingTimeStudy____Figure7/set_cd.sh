vec=(-1.55
-1.16
-0.8
-0.4
-0.2
-0.1
-0.06
0.0
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
1.55
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
