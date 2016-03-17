#!/usr/bin/env bash

rootdir=~/data/gfp_pirna/data/downloaded
for i in Zamore.RSQ.fly.ago3_mut_aub_mut.rep1.SRR1924207   Zamore.RSQ.fly.ago3_mut_piwi_mut.rep2.SRR1924210  Zamore.RSQ.fly.aub_mut.rep1.SRR1924203   Zamore.RSQ.fly.piwi_mut.rep2.SRR1924202 Zamore.RSQ.fly.ago3_mut_aub_mut.rep2.SRR1924208   Zamore.RSQ.fly.ago3_mut.rep1.SRR1924205           Zamore.RSQ.fly.aub_mut.rep2.SRR1924204   Zamore.RSQ.fly.WT.rep1.SRR1783725 Zamore.RSQ.fly.ago3_mut_piwi_mut.rep1.SRR1924209  Zamore.RSQ.fly.ago3_mut.rep2.SRR1924206           Zamore.RSQ.fly.piwi_mut.rep1.SRR1924201  Zamore.RSQ.fly.WT.rep2.SRR1783726; do    
    nf=$(grep 'unique_mapper_reads' ${rootdir}/${i}/*.basic_stats | cut -d $'\t' -f2 | awk '{ print 1e6 / $1 }')
    bam=$(readlink -e ${rootdir}/${i}/genome_mapping/*.sorted.bam)
    ## echo ${bam}
    echo "bedtools bamtobed -bed12 -tag NH -i ${bam} | awk -v strand=1 'BEGIN{FS=OFS=\"\t\"}{e=substr(\$4,length(\$4)); if (e==strand) \$6=(\$6==\"+\"?\"-\":\"+\"); if (\$5==1) print \$0; }' > ${i}.sorted.unique.bed12 "
    echo "bash count_transposon_groups.sh ${i}.sorted.unique.bed12 ${nf}"
    # if [[ ! -f ${bam} ]]; then
    # 	echo "${bam} is not ready"
    # fi
    # run1.sh ${rootdir}/${i}/genome_mapping/ ${nf}
done
