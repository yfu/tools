#!/usr/bin/env python

# Get some basic stats for your gff file
# Notice that the upstream1000, upstream2000 and upstream5000 files
# may cross chromosome boundaries

import sys
import gffutils
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# fn = "small.gff"
# fn = "hi5_v0.8.2.run6.gene.gff"
# fasta_fn = "hi5_v0.8.2.fasta"
fasta_fn = sys.argv[1]
fn = sys.argv[2]
fasta = Fasta(fasta_fn)

# output="test"
output = sys.argv[3]

output_cds_fn          = output + ".cds.fa"
output_pep_fn          = output + ".pep.fa"
output_mrna_fn         = output + ".mRNA.fa"
output_upstream1000_fn = output + ".upstream1000.bed"
output_upstream2000_fn = output + ".upstream2000.bed"
output_upstream5000_fn = output + ".upstream5000.bed"

output_cds          = open(output_cds_fn, "w")
output_pep          = open(output_pep_fn, "w")
output_mrna         = open(output_mrna_fn, "w")
output_upstream1000 = open(output_upstream1000_fn, "w")
output_upstream2000 = open(output_upstream2000_fn, "w")
output_upstream5000 = open(output_upstream5000_fn, "w")

# fn = gffutils.example_filename('intro_docs_example.gff')
db = gffutils.create_db(fn, dbfn=fn + ".db",
                        force=True, keep_order=True,
                        merge_strategy='merge',
                        sort_attribute_values=True)



for mrna in db.features_of_type("mRNA", order_by="start"):
    print mrna
    merged_cds = ""
    for cds in db.children(mrna, featuretype="CDS", order_by='start'):
        print cds
        merged_cds += str(cds.sequence(fasta, use_strand=True))
        strand = cds.strand
    print merged_cds
    merged_cds = Seq(merged_cds, generic_dna)

    merged_exon = ""
    for exon in db.children(mrna, featuretype="exon", order_by="start"):
        print exon
        merged_exon += str(exon.sequence(fasta, use_strand=True))
        strand = exon.strand
    print merged_exon
    merged_exon = Seq(merged_exon, generic_dna)

    if strand == "-":
        merged_cds = merged_cds.reverse_complement()
        merged_exon = merged_exon.reverse_complement()
    pep = merged_cds.translate()

    assert len(mrna["Name"]) == 1
    print >>output_pep, ">" + mrna["Name"][0]
    print >>output_pep, pep

    print >>output_cds, ">" + mrna["Name"][0]
    print >>output_cds, merged_cds

    print >>output_mrna, ">" + mrna["Name"][0]
    print >>output_mrna, merged_exon

for gene in db.all_features(featuretype="gene", order_by="start"):
    assert len( list(gene["Name"]) ) == 1
    gene_name = gene["Name"][0]
    start = gene.start - 1
    end = gene.end
    strand = gene.strand
    seqid = gene.seqid
    if strand == "+":
        print >>output_upstream1000, "\t".join(
            (seqid, str(start-1000), str(start),
             gene_name + "_up1000", ".", "+")
        )
        print >>output_upstream2000, "\t".join(
            (seqid, str(start-2000), str(start),
             gene_name + "_up2000", ".", "+")
        )
        print >>output_upstream5000, "\t".join(
            (seqid, str(start-5000), str(start),
             gene_name + "_up5000", ".", "+")
        )
    elif strand == "-":
        print >>output_upstream1000, "\t".join(
            (seqid, str(end), str(end + 1000),
             gene_name + "_up1000", ".", "-")
        )
        print >>output_upstream2000, "\t".join(
            (seqid, str(end), str(end + 2000),
             gene_name + "_up2000", ".", "-")
        )
        print >>output_upstream5000, "\t".join(
            (seqid, str(end), str(end + 5000),
             gene_name + "_up5000", ".", "-")
        )

