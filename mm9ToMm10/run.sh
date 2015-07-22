
for i in prepachytene hybrid pachytene; do
    bedToGenePred piRNA.cluster.${i}.bed12 stdout | genePredToGtf file stdin stdout | sed -E 's/gene_id \"([^;]+)\.[0-9]"/gene_id \"\1\"/g' | sed -E 's/gene_id ([^ ]+)/gene_id \1 gene_name \1/' > piRNA.cluster.${i}.gtf
done
