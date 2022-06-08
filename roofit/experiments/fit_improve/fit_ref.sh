#!/bin/bash
# @Author: Billy Li <billyli>
# @Date:   05-19-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 05-19-2022

fullFileDir="../../../analysis/allEvents"
for i in $fullFileDir/fullWR*.root; do
  filename=$(basename $i)
  root -l -b RooExpmCB.cxx "ref.C(\"$filename\")" -q
done
