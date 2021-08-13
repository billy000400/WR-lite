# @Author: Billy Li <billyli>
# @Date:   08-05-2021
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 08-12-2021



addressTwo="/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/"
address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR800_N400/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR800_N400/out_WR800N400_${i}.root isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR800_N600/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR800_N600/out_WR800N600_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR800_N700/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR800_N700/out_WR800N700_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


#Start WR 1000

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1000_N200/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1000_N200/out_WR1000N200_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1000_N900/"
for i in {1..28}
do
	if (($i!=19)); then
		cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1000_N900/out_WR1000N900_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
		if (( $i%7==0 )); then
			wait
		fi
	fi
done


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1000_N400/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1000_N400/out_WR1000N400_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1000_N600/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1000_N600/out_WR1000N600_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1000_N800/"
for i in {1..28}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1000_N800/out_WR1000N800_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done



address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR800_N200/"
for i in {1..6}
do
	cmsRun python/ExtractRecoMass_WR_N.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR800_N200/out_WR800N200_${i}.root  isSignal=True genTrainData=True  trainFile=ml.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done
