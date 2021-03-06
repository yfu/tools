#!/usr/bin/env bash

for D in SRS.GS1CT1.br1.tr1.ox.ovary SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do
    # for D in SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do
    # D=SRS.GS1CT1.br1.tr1.ox.ovary
    PREFIX=$D
    GCR_D=${D}/gene+cluster+repBase
    GCR_ALLMAP_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    PREFIX_COUNT=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    
    echo "awk '{ if(\$1 ~ /cluster|42AB|flam/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.cluster.count"
    
    echo "awk '{ if(\$1 ~ /^FBgn/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.repBase.count"
    
    echo "awk '{ if(\$1~/^NM_|^NR_/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.gene.count"
    
    echo "awk '{if(\$1~/nos-gal4-vp16|GSV6|UASp-EGFP/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.construct.count"
done
