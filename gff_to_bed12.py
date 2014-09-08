#!/usr/bin/env python
"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# This script will convert a GFF file into a bed12 file, which describes the exon information

# Author: Yu Fu

import sys, re

class Exon:
    """A simple class representing one exon in a GFF file"""
    def __init__(self, chrom, start, end, strand, gene):
        self.chrom, self.start, self.end, self.strand, self.gene = chrom, start, end, strand, gene
    def __repr__(self):
        return "Exon(" + self.chrom + ":" + str(self.start) + "-" + str(self.end) + "\t" + self.strand + "\t" + self.gene + ")"

def get_bed12(exons):
    """Given a set of exons (belonging to the same gene), print one bed12 entry"""
    """ Example
chr1    11873   14409   uc001aaa.3      0       +       11873   11873   0       3       354,109,1189,   0,739,1347,
chr1    11873   14409   uc010nxq.1      0       +       12189   13639   0       3       354,127,1007,   0,721,1529,
chr1    11873   14409   uc010nxr.1      0       +       11873   11873   0       3       354,52,1189,    0,772,1347,
"""
    start = sys.maxint
    end = 0
    abs_starts = []
    abs_ends = []
    for e in exons:
        chrom = e.chrom
        abs_starts.append(e.start)
        abs_ends.append(e.end)
        if e.start < start:
            start = e.start
        if e.end > end:
            end = e.end
        block_count = len(abs_starts)
        block_sizes = ""
        block_starts = ""
        for i in range(len(abs_starts)):
            block_sizes += str(abs_ends[i] - abs_starts[i]) + ","
            block_starts += str(abs_starts[i] - start) + ","
    return "\t" . join([ chrom, str(start), str(end), e.gene, "0", e.strand, str(start), str(start), "0", str(block_count), block_sizes, block_starts ])
    
def main():
    if len(sys.argv) <= 1:
        print "Usage: gff_to_bed12.py test_data/TAIR10_GFF3_genes_test.gff"
        sys.exit(0)
    # Example: Parent=AT1G01010.1
    pat = re.compile(r"Parent=(.+)")
    fh = open(sys.argv[1])

    # The value is a list, containing one or more exons
    genes = {}

    for line in fh.readlines():
        line = line.strip()
        e = line.split()
        [chrom, source, feature, start, end, score, strand, frame, attribute] = e
        start = int(start) - 1      # Notice that bed is 1-based and gff is 0-based
        # Notice that in a bed file, end position is not
        # included, whereas in a gff file, the end is included
        end = int(end)

        if feature != "exon":
            continue
        mat = pat.search(attribute)
        gene_name = mat.group(1)
        assert gene_name is not None

        # Just in case we encounter bad GFF file
        # print Exon(chr, start, end, strand, gene_name)
        if gene_name not in genes:
            genes[gene_name] = [ Exon(chrom, start, end, strand, gene_name) ]
        else:
            genes[gene_name].append( Exon(chrom, start, end, strand, gene_name) )

    fh.close()
    for g in genes:
        # print genes[g]
        exons = genes[g]
        print get_bed12(exons)

if __name__ == "__main__":
    main()
