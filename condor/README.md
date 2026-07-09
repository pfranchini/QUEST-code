
Some scripts for the HTCondor batch system (Oxford/RHUL, mostly used now in Oxford)

`batch`: contains scripts used for Oxford ND3 simulation, both for non-rethrow and rethrow volumes.

Usage:
```
   ./submit.sh <start> <stop>
```

Merge root files (in vanilla shell):
```
   hadd -dbg -f -k merge.root */*.root
```   
