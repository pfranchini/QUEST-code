dir=$1

if [ ! -d "$dir" ]; then
   exit 1
fi
 
cd $dir/output

echo "Counting and merging..."
roots=`find . -type f -name "output*.root" -size +1k -exec ls -lh {} \; | wc -l`
jobs=`grep -R "End of the session" * | wc -l`

hadd -v 0 -f -k -j 10 merge.root */*.root
events=`grep -h "Number of events" 1/* | tail -n1 | awk {'print $4'}`

echo
echo "Root files completed: " $roots
echo "Jobs completed:       " $jobs
echo "Simulated events/job: " $events

#root -q ~/dataQUEST/QUEST/ND3/batch/code/plot.c



cd -
