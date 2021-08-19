# @Author: Billy Li <billyli>
# @Date:   08-05-2021
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 08-17-2021



addressTwo="/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/"

#Start WR 2000

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR6000_N3000/"
for i in {1..19}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR6000_N3000/out_WR6000N3000_${i}.root isSignal=True   > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done
