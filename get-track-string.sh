# /data/www/zlab-trackhub/yfu/human_pachytene/best_aln

# Script name
sn=$(basename "$0")
if [[ "$#" -eq 0 ]]; then
    echo "Usage: $sn your.awesome.bigWig" 1>&2
    echo "       $sn your.awesome.bigBed" 1>&2
    echo "       $sn your.awesome.Watson.bigWig" 1>&2
    echo "Do not use *.Crick.bigWig" 1>&2
    exit 1
fi

url_root=zlab-trackhub.umassmed.edu
p=$1
p=$(readlink -e $p)
rel=${p#/data/www/zlab-trackhub/}
rel=$(dirname ${rel})
echo ${rel}
fn=$(basename ${p})
suf=${fn##*.}
if [[ $suf == "bigwig" || $suf == "bigWig" ]]; then
    type=bigWig
elif [[ $suf == "bigbed" || $suf == "bigBed" ]]; then
    type=bigBed
else
    echo "Unknown file extension: $suf"
fi
    
prefix=${i%.Watson.bigWig}

if [[ $fn =~ ".*Watson.*" ]]; then
    # Signals with strand
    crick=$(cat ${fn} | sed 's/Watson/Crick/')
    echo "track type=${type} name=\"${fn}\" bigDataUrl=http://${url_root}/${rel}/${fn} color=0,0,255 visibility=2 priority=50"
    echo "track type=${type} name=\"${crick}\" bigDataUrl=http://${url_root}/${rel}/${crick} color=255,0,0 visibility=2 priority=50"
else
    # Signals without strand
    echo "track type=${type} name=\"${fn}\" bigDataUrl=http://${url_root}/${rel}/${fn} color=0,0,255 visibility=2 priority=50"
fi





