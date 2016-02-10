#!/usr/bin/env python

# I have given up this method. It cannot extract all info (query_id, or the name of the query sequence) from the XML file
from Bio.Blast import NCBIXML

result_handle = open("SRA_to_DEG_flanking.test100.xml")
# result_handle = open("SRA_to_DEG_flanking.xml")
blast_records = NCBIXML.parse(result_handle)
# E_VALUE_THRESH = 0.04
E_VALUE_THRESH = 1000
DEBUG=False



# print  hsp.query, hsp.sbjct, hsp.match, query_name, sbjct_name, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.align_length, hsp.expect, hsp.score, strand
print "# query_match_seq\tsbjct_match_seq\tmatch\tquery_name\tsbjct_name\thsp.query_start\thsp.query_end\thsp.sbjct_start\thsp.sbjct_end\talign_length\teval\tbitscore\tstrand"

for blast_record in blast_records:
    if DEBUG == True:
        print "-" * 80
    # print dir(blast_record)
    # This is missing in the document...
    query_name = blast_record.query

    tmp = query_name.split("_")
    query_copy = tmp[0]
    query_seq = tmp[1]
    # print "blast_record\t" blast_record.query
    for alignment in blast_record.alignments:
        # print dir(alignment)
        # print vars(alignment)
        sbjct_name = alignment.hit_def
        # print alignment.query
        # print('****Alignment****')
        # print 'sequence', alignment
        # print('length:', alignment.length)
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                # print vars(hsp)
                # print dir(hsp)
                # print('positives', hsp.positives)
                # print('identities', hsp.positives)
                # print('e value:', hsp.expect)
                # print('score:', hsp.score)
                # print(hsp.query[0:75] + '...')
                # print(hsp.match[0:75] + '...')
                # print(hsp.sbjct[0:75] + '...')
                # print hsp.query
                # print hsp.match
                # print hsp.sbjct
                if hsp.sbjct_end < hsp.sbjct_start:
                    strand = "-"
                else:
                    strand = "+"
                # Use sbjct as the reference
                n_mat = 0
                n_mis = 0
                n_ins = 0
                n_del = 0
                for i in range(0, len(hsp.query)):
                    if hsp.query[i] == hsp.sbjct[i]:
                        n_mat += 1
                    elif hsp.query[i] == "-":
                        n_ins += 1
                    elif hsp.sbjct[i] == "-":
                        n_del += 1
                    else:
                        n_mis += 1
                if DEBUG == True:
                    print "# matches\t" + str(n_mat)
                    print "# mismatches\t" + str(n_mis)
                    print "# insertions\t" + str(n_ins)
                    print "# deletions\t" + str(n_del)
                    print "strand\t" + strand
                # print(hsp.sbjct)
                # hsp.strand is not available for blastn XML files
                # print hsp.strand
                # print dir(hsp)
                print "\t".join( (hsp.query, hsp.sbjct, hsp.match, query_name, sbjct_name, str(hsp.query_start), str(hsp.query_end), str(hsp.sbjct_start), str(hsp.sbjct_end), str(hsp.align_length), str(hsp.expect), str(hsp.score), strand) )
