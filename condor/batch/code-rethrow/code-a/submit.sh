# Submit stop-start batches of a number=Queue of identical jobs as specified in job.submit
#
# Usage: ./submit.sh <start> <stop> <volume> <ion>
#

start=$1
stop=$2
volume=$3
ion=$4

dir=/data/questdmc/users/franchinip/QUEST/ND3/${volume}-${ion}/output-a
mkdir -p $dir

echo $dir

read -p "Are you sure? [y/n] " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

    echo "Submitting jobs..."
    
    for ((i = $start ; i <= $stop ; i++)); do
	echo "condor_submit job.submit job=$i"
	condor_submit job.submit job=$i volume=$volume ion=$ion
	sleep 0.5
    done
    
    sleep 2
    condor_q

fi
