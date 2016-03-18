#!/usr/bin/env bash

# bash run4-collect-z1.sh > all.z1
for i in SRS.GS1BT1.br1.tr1.ox.ovary  SRS.GS1BT1.br2.tr1.ox.ovary  SRS.GS1BT1.br3.tr1.ox.ovary; do
    for j in $(seq 0 500 10000); do
	## SRS.GS1BT1.br2.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0.uniq.piRNA.uniq.GSV6_downstream_in_zip.dm3.2k.offset4000.bed.Watson3to5.phasing_z1
	offset=-${j}
	up=${i}.x_rRNA.x_hairpin.dm3_chengjianv0.uniq.piRNA.uniq.GSV6_upstream_in_zip.dm3.2k.offset${offset}.bed
	z1=$(cat ${up}.Watson3to5.phasing_z1)
	echo -e "${i}\t${offset}\tWatson\t${z1}"
	z1=$(cat ${up}.Crick5to3.phasing_z1)
	echo -e "${i}\t${offset}\tCrick\t${z1}"
	
	offset=${j}
	down=${i}.x_rRNA.x_hairpin.dm3_chengjianv0.uniq.piRNA.uniq.GSV6_downstream_in_zip.dm3.2k.offset${offset}.bed
	z1=$(cat ${down}.Watson3to5.phasing_z1)
	echo -e "${i}\t${offset}\tWatson\t${z1}"
	z1=$(cat ${down}.Crick5to3.phasing_z1)
	echo -e "${i}\t${offset}\tCrick\t${z1}"
    done
done



