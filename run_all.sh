#!/usr/bin/env bash

# Just a script to submit lots of jobs
for i in HEIL JMAC KNBD OQSU PRTV; do
    for j in id1 id2 id3 id4 id5 id6; do
	cat << "EOF" >> test_qsub.txt
#!/usr/bin/env bash
#$ -cwd
#$ -o $HOME/log/$JOB_ID.out
#$ -e $HOME/log/$JOB_ID.err
#$ -M fuyu12345@gmail.com
#$ -S /bin/bash
#$ -V
###$ -pe single 4
#$ -l mem_free=15G
###$ -q infiniband.q
EOF
# The above part does not need parameter substitution
# whereas the following part does
cat <<EOF
cd ${i}_${j}
~hanb/nearline/PE_Pipeline/bin/run.sh -l ${i}.r1_${j}.fq -r ${i}.r2_${j}.fq -s fly -c
EOF


    done;
done;
