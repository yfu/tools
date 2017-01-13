
if [[ "$#" -ne 1 ]]; then
    echo "Given a maf file and a speices name, this script reports the total alignment length. "
    echo "You probably want to filter the maf file first as some maf files contain single-row alignment"
    echo "Usage: mafAlnLen.sh your_species_name < input.maf"
    exit 1
fi

sname=$1
awk -v FS=" " -v OFS=" " -v sname=dmel '{ split($2, a, "."); if(a[1]==sname) { s+=$4} } END{ print s }' /dev/stdin
