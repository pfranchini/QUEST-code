dir=$1

#=============================
check_existing=true
#=============================

if [ ! -d "$dir" ]; then
   exit 1
fi
 
cd $dir/output-b

echo "Counting..."
roots=`find . -type f -name "output*.root" -size +1k -exec ls -lh {} \; | wc -l`
jobs=$(find . -type f -name '*.out' -exec grep -H "End of the session" {} + | wc -l)

merge=true
if [ "$check_existing" = true ] ; then
    if [ -f "merge.root" ]; then
	filesize=$(stat -c%s "merge.root")
	if [ "$filesize" -gt 10485 ]; then
	    echo "merge.root is already present."
	    ls -l merge.root
	    merge=false
	fi
    fi
fi
	
if [ "$merge" = true ] ; then
    echo "Merge files..."
    #hadd -v 0 -f -k -j 10 merge.root */*.root
    hadd -v 0 -f -k -j 5 merge.root */*.root
fi

#hadd -v 0 -f -k merge.root 1/*.root  # otherwise the output is too big
events=`grep -h "Number of events" 1/*.out | tail -n1 | awk {'print $4'}`

echo
echo "Root files completed: " $roots
echo "Jobs completed:       " $jobs
echo "Simulated events/job: " $events

#root -q ~/dataQUEST/QUEST/ND3/batch/code-rethrow/code-a/plot-rethrow.c

cd -
