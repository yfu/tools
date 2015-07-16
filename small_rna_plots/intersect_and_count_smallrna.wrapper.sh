#!/usr/bin/env bash



cmd=parallel.commands
if [[ -f ${cmd} ]]; then
    rm ${cmd}
fi

rootdir=/home/fuy2/data/gfp_pirna/data
for i in SRS.42A18AT1.br1.tr1.ox.ovary.fq.gz SRS.42A18AT1.br2.tr1.ox.ovary.fq.gz SRS.42A18AT1.br3.tr1.ox.ovary.fq.gz SRS.GS1AT1.br1.tr1.ox.ovary.fq.gz SRS.GS1AT1.br2.tr1.ox.ovary.fq.gz SRS.GS1AT1.br3.tr1.ox.ovary.fq.gz SRS.42A18BT1.br1.tr1_2.ox.ovary.fq.gz SRS.42A18BT1.br2.tr1.ox.ovary.fq.gz SRS.42A18BT1.br3.tr1.ox.ovary.fq.gz SRS.GS1BT1.br1.tr1.ox.ovary.fq.gz SRS.GS1BT1.br2.tr1.ox.ovary.fq.gz SRS.GS1BT1.br3.tr1.ox.ovary.fq.gz SRS.42A18CT1.br1.tr1.ox.ovary.fq.gz SRS.42A18CT1.br2.tr1.ox.ovary.fq.gz SRS.42A18CT1.br3.tr1.ox.ovary.fq.gz SRS.GS1CT1.br1.tr1.ox.ovary.fq.gz SRS.GS1CT1.br2.tr1.ox.ovary.fq.gz SRS.GS1CT1.br3.tr1_2.ox.ovary.fq.gz; do
    prefix=${i%.fq.gz}
    bash intersect_and_count_smallrna.sh -N allxmirna -a ${rootdir}/${prefix} -o ${cmd} -g dm3_chengjian -A ${prefix}
done
