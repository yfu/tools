for i in prepachytene hybrid pachytene; do
    mkdir -p ${i}
done

for i in prepachytene hybrid pachytene; do
    cluster_mm10=~/data/shared/mm9ToMm10/piRNA.cluster.${i}.bed6
    prefix_mm10=$(basename ${cluster_mm10%.bed12})
    prefix_mm10=${prefix_mm10%.bed}.mm10
    cluster_watson_mm10=${prefix_mm10}.Watson.bed6
    cluster_crick_mm10=${prefix_mm10}.Crick.bed6
    cluster_nostrand_mm10=${prefix_mm10}.nostrand.bed6
    # Make the name of each range unique, required by bigWigAverageOverBed
    awk 'BEGIN{FS=OFS="\t"; n=1} { if($6=="+") {$4=$4 "_" n; print $0; n+=1} }' ${cluster_mm10} > ${cluster_watson_mm10}
    awk 'BEGIN{FS=OFS="\t"; n=1} { if($6=="-") {$4=$4 "_" n; print $0; n+=1} }' ${cluster_mm10} > ${cluster_crick_mm10}
    awk 'BEGIN{FS=OFS="\t"; n=1} { $4=$4 "_" n; print $0; n+=1 }' ${cluster_mm10} > ${cluster_nostrand_mm10}

    cluster_mm9=~/data/shared/mm9/piRNA.cluster.${i}.bed6
    prefix_mm9=$(basename ${cluster_mm9%.bed12})
    prefix_mm9=${prefix_mm9%.bed}.mm9
    cluster_watson_mm9=${prefix_mm9}.Watson.bed6
    cluster_crick_mm9=${prefix_mm9}.Crick.bed6
    cluster_nostrand_mm9=${prefix_mm9}.nostrand.bed6
    # Make the name of each range unique, required by bigWigAverageOverBed
    awk 'BEGIN{FS=OFS="\t"; n=1} { if($6=="+") {$4=$4 "_" n; print $0; n+=1} }' ${cluster_mm9} > ${cluster_watson_mm9}
    awk 'BEGIN{FS=OFS="\t"; n=1} { if($6=="-") {$4=$4 "_" n; print $0; n+=1} }' ${cluster_mm9} > ${cluster_crick_mm9}
    awk 'BEGIN{FS=OFS="\t"; n=1} { $4=$4 "_" n; print $0; n+=1 }' ${cluster_mm9} > ${cluster_nostrand_mm9}

    for j in $(ls -1 ../../data/*signal*bigWig); do
    	bigwig=${j}
    	# echo ${prefix}
	
    	track=$(basename ${j%.bigWig})
    	signal_watson=${i}/${track}.${i}_exons.Watson.signal
    	signal_crick=${i}/${track}.${i}_exons.Crick.signal
    	signal=${i}/${track}.${i}.nostrand.signal # For the case where the bigWig file does not specify the strand
	if [[ $j =~ mm10 ]]; then
    	    if [[ $j =~ plus ]]; then
    		# If the signals have both Watson and Crick strand, average the separately
    		echo bigWigAverageOverBed ${bigwig} ${cluster_watson_mm10} ${signal_watson}
    	    elif [[ $j =~ minus ]]; then
    		echo bigWigAverageOverBed ${bigwig} ${cluster_crick_mm10} ${signal_crick}
    	    else	
    		echo bigWigAverageOverBed ${bigwig} ${cluster_nostrand_mm10} ${signal}
    	    fi
	elif [[ $j =~ mm9 ]]; then
    	    if [[ $j =~ plus ]]; then
    		# If the signals have both Watson and Crick strand, average the separately
    		echo bigWigAverageOverBed ${bigwig} ${cluster_watson_mm9} ${signal_watson}
    	    elif [[ $j =~ minus ]]; then
    		echo bigWigAverageOverBed ${bigwig} ${cluster_crick_mm9} ${signal_crick}
    	    else	
    		echo bigWigAverageOverBed ${bigwig} ${cluster_nostrand_mm9} ${signal}
    	    fi
	else
	    echo "Warning! Found one bigWig file without genome assembly specified" >/dev/stderr
	fi
    done
done
