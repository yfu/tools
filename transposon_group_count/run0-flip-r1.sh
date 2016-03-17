#!/usr/bin/env bash


root=~/data/gfp_pirna/data/
# for i in SRS.AC1.br1.tr1.ox.ovary SRS.42A18AT1.br1.tr1.ox.ovary SRS.42A18AT1.br2.tr1.ox.ovary SRS.42A18AT1.br3.tr1.ox.ovary SRS.42A18AT1.br4.tr1_2.ox.ovary SRS.GS1AT1.br1.tr1.ox.ovary SRS.GS1AT1.br1.tr1.unox.ovary SRS.GS1AT1.br2.tr1.ox.ovary SRS.GS1AT1.br3.tr1.ox.ovary SRS.ZIPAT1.br1.tr1.ox.ovary ; do
for i in     RSRZ.CC1.br{1,2,3}.tr1.ovary RSRZ.GS1CT1.br{1,2,3}.tr1.ovary \
			RSRZ.CC2.br{1,2,3}.tr1.ovary RSRZ.GS1CT2.br{1,2,3}.tr1.ovary \
			RSRZ.CC3.br{1,2,3}.tr1.ovary RSRZ.GS1CT3.br{1,2,3}.tr1.ovary \
			RSRZ.CC4.br{1,2,3}.tr1.ovary RSRZ.GS1CT4.br{1,2,3}.tr1.ovary \
			RSRZ.AC1.br0.tr1.ovary \
			RSRZ.GS1AT1.br0.tr1.ovary \
			RSRZ.42A18AT1.br0.tr1.ovary \
			RSRZ.ZIPAT1.br0.tr1.ovary; do
    bam=$(readlink -e ${root}/${i}/genome_mapping/*.sorted.bam)
    echo "bam_to_bed12_flip_read1.sh ${bam} > ${i}.bed12"
done
    
    
