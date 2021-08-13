# @Author: Billy Li <billyli>
# @Date:   08-05-2021
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 08-12-2021



addressTwo="/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/"

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1400_N1000/"
for i in {1..28}
do
	if (($i!=2)); then
		cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1400_N1000/out_WR1400N1000_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
		if (( $i%7==0 )); then
			wait
		fi
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1400_N1200/"
for i in {1..27}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1400_N1200/out_WR1400N1200_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


#Start WR1600

wait


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N400/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N400/out_WR1600N400_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N1500/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N1500/out_WR1600N1500_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N600/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N600/out_WR1600N600_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done


address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N800/"
for i in {1..28}
do
	if (($i!=11)); then
		cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N800/out_WR1600N800_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
		if (( $i%7==0 )); then
			wait
		fi
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N1000/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N1000/out_WR1600N1000_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N1200/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N1200/out_WR1600N1200_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done

address="/hdfs/cms/user/krohn045/WR_SignalSamples/WR1600_N1400/"
for i in {1..28}
do
	cmsRun python/cfg.py inputFiles=file:${address}MINIAOD_${i}.root  outputFile=${addressTwo}WR1600_N1400/out_WR1600N1400_${i}.root isSignal=True genTrainData=True  trainFile=ml3.txt > /dev/null &
	if (( $i%7==0 )); then
		wait
	fi
done
