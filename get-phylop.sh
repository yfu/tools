if [[ "$#" == 0 ]]; then
    echo "Usage: get-phylop.sh your.bed"
    exit
fi
cpu=32
fn=$1
i=1
d=tmp_query_phylop.${RANDOM}${RANDOM}
mkdir -p ${d}
commands=${RANDOM}${RANDOM}.commands
while read line; do
    echo "${line}" > ${d}/${i}.bed
    echo "bed_to_phylop_one_line.sh ${d}/${i}.bed"
    i=$(( i + 1 ))
done < ${fn} > ${commands}
cat ${commands} | parallel --progress -j ${cpu}
ls ${d}/*.with_phylop | xargs cat | sort -k1,1 -k2,2n - > ${fn}.phylop
rm -r ${d}
