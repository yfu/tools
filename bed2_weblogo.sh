

input=$1
cat ${input} | awk '{ d[$7]+=$4/$5 } END{ for(i in d) { print i "\t" d[i] } }' | awk '{ for(i=0; i<$2; i++) { print ">"++j; print substr($1, 1, 23) } }' | weblogo -F pdf -P "" > ${input}.weblogo.entropy.pdf
cat ${input} | awk '{ d[$7]+=$4/$5 } END{ for(i in d) { print i "\t" d[i] } }' | awk '{ for(i=0; i<$2; i++) { print ">"++j; print substr($1, 1, 23) } }' | weblogo -F pdf -P "" -U probability > ${input}.weblogo.probability.pdf
