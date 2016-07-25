### BSUB -q short # which queue we want to run in  

submit_jobs() {
    while read line; do
	cmd=${line}
	echo "
#!/bin/bash

#BSUB -n 2
#BSUB -R rusage[mem=50000] # ask for 2GB per job slot, or 8GB total
#BSUB -W 2:00
#BSUB -q short # which queue we want to run in
#BSUB -R \"span[hosts=1]\" # All job slots on the same node (needed for threaded applications)
#BSUB -e logs/out.%J.%I.err
#BSUB -o logs/out.%J.%I.out
export PATH=$HOME/dazzler/DALIGNER:$HOME/dazzler/DAZZ_DB:$HOME/dazzler/DEXTRACTOR:$PATH:$PATH
$cmd
" | bsub
 
    done < $1
}
mkdir -p logs

submit_jobs jobs.sh




