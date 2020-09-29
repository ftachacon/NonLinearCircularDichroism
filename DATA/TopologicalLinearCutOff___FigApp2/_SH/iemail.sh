
cd ./$3 

mail -n -s "phi = $1 - $2 inter" alexis.chacon@lanl.gov < interband_dipole_full_evol.dat 
mail -n -s "phi = $1 - $2 intra" alexis.chacon@lanl.gov < intraband_current_full_evol.dat 
mail -n -s "phi = $1 - $2 laser" alexis.chacon@lanl.gov < outlaserdata.dat 
mail -n -s "phi = $1 - $2 params" alexis.chacon@lanl.gov < setOfparameters.dat 

pwd 

echo ""
echo "phi= $1;   ellip = $2"
echo ""

cd ../
