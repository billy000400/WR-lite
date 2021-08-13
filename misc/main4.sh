# @Author: Billy Li <billyli>
# @Date:   08-05-2021
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 08-12-2021



addressTwo="/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/"
#Start WR 1800

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N400/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N400/out_WR1800N400_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N1700/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N1700/out_WR1800N1700_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N600/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N600/out_WR1800N600_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N800/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N800/out_WR1800N800_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N1000/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N1000/out_WR1800N1000_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N1200/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N1200/out_WR1800N1200_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N1400/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N1400/out_WR1800N1400_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1800_N1600/"
for i in {1..27}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1800_N1600/out_WR1800N1600_${i}.root isSignal=True genTrainData=True  trainFile=ml4.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done
