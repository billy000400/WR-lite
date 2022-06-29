# @Author: Billy Li <billyli>
# @Date:   06-29-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 06-29-2022

sampleFileDir="../../../data/RooFitMC/ratio_1e-2/WR1600N800/"
declare -i i=1
for f in $sampleFileDir/*.root; do
  filename=$(basename $f)
  root -l -b -q RooExpm.cxx RooExpmCB.cxx "fitMassRatio.C(\"$filename\",\"1e-2\")" > /dev/null &
  echo "processing $filename"
  if (( $i%8==0 )); then
      wait
  fi
  i=$(($i + 1))
done
