
# Usage: ./submit.sh <start> <stop>

start=$1
stop=$2

cat job.submit | grep basedir | grep scratch

read -p "Are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

    echo "Submitting jobs..."
    
    for ((i = $start ; i <= $stop ; i++)); do
	echo "condor_submit job.submit job=$i"
	condor_submit job.submit job=$i
	sleep 0.5
    done
    
    sleep 2
    condor_q

fi
