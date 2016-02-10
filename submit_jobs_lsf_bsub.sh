
# awk '{ if($0 ~ /^#/) { f+=1 } else { print $0 > "run_daligner.part" f ".sh" } }' run_daligner.sh

### BSUB -q short # which queue we want to run in  

submit_jobs() {
    while read line; do
        cmd=${line}
        echo "
#!/bin/bash

#BSUB -n 4
#BSUB -R rusage[mem=64000]
#BSUB -W 168:00
#BSUB -q long # which queue we want to run in
#BSUB -R \"span[hosts=1]\" # All job slots on the same node (needed for threaded applications)
#BSUB -e logs/out.%J.%I.err
#BSUB -o logs/out.%J.%I.out
export PATH=$HOME/dazzler/DALIGNER:$HOME/dazzler/DAZZ_DB:$HOME/dazzler/DEXTRACTOR:$PATH:$PATH
$cmd
" | bsub
 
    done < $1
}
mkdir -p logs
# submit_jobs run_daligner_preads.part1.sh
# submit_jobs run_daligner_preads.part2.sh
# submit_jobs cmds.sh
submit_jobs cmds_graph_to_contig.sh
## submit_jobs test_run_daligner_preads.part1.sh
