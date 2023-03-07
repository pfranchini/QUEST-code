#!/bin/sh

#==========================================

process=$1
job=$2
dir=$3

code=/home/pfranchi/scratch/QUEST/cosmics/fridge4/QUEST-detector-simulation/build
events=100 #1000000

#==========================================

output_file=${dir}/${job}/output-${process}.root

cd ${code}
source ${code}/env.sh
#${code}/sim -p cosmic -n ${events} --stepping -V -z --output ${output_file}
${code}/g4quest -p cosmic -n ${events} -z --output ${output_file}


