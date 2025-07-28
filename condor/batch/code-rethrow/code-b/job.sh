#!/bin/sh

#==========================================

process=$1
job=$2
dir=$3  #/data/questdmc/users/franchinip/QUEST/ND3/$(volume)-$(ion)/output-b
volume=$4
ion=$5

code=/data/questdmc/users/franchinip/QUEST/ND3/QUEST-detector-simulation/build
events=1000000 #1M

#==========================================

output_file=${dir}/${job}/output-${process}.root

cd ${code}
#source ${code}/env.sh
export QT_QPA_PLATFORM=offscreen

echo "Process:" $process

# to avoid duplicate Geant seeds
sleep $((process * 2))

# g4quest command line:
${code}/g4quest -f ND3 -z --generator rethrow --input /data/questdmc/users/franchinip/QUEST/ND3/${volume}-${ion}/output-a/merge.root -n ${events} --output ${output_file}

