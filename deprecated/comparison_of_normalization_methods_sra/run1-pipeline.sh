#!/usr/bin/env bash

data_root=../../data/
for i in SRS.GS1CT1.br1.tr1.ox.ovary.fq.gz SRS.GS1CT1.br2.tr1.ox.ovary.fq.gz SRS.GS1CT1.br3.tr1_2.ox.ovary.fq.gz SRS.GS1CT1.br3.tr1.ox.ovary.fq.gz SRS.GS1CT1.br3.tr2.ox.ovary.fq.gz; do
    prefix=${i%.fq.gz}
    log=${prefix}.log
    output=${prefix}
    echo "piPipes small -i ${data_root}/${i} -g dm3_chengjian -o ${output} -c 16 &> ${log}"
done
