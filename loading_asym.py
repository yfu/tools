#!/usr/bin/env python3

## Given a bed2 file, this script checkes the loading asymmetry of siRNAs

import sys
import argparse

def loading_asym(fp_watson, fp_crick, tp_watson, tp_crick):
    # Pseudocount to avoid dividing by 0
    pc = 0.1
    print("# chrom\tthree_prime_pos_on_watson\tasym_score")
    for chrom in tp_watson:
        for end in tp_watson[chrom]:

            if end - 2 in fp_crick[chrom]:
                partner = fp_crick[chrom][end-2]
            else:
                partner = 0
            asym_score = (tp_watson[chrom][end] + pc) / (partner + pc)
            print("{}\t{}\t{}".format(chrom, end, asym_score))

            
def loading_asym_rev(fp_watson, fp_crick, tp_watson, tp_crick):
    # Just for testing
    # If I traverse the crick strand, I should get the same results
    # Pseudocount to avoid dividing by 0
    pc = 0.1
    print("# chrom\tthree_prime_pos_on_watson\tasym_score")
    for chrom in tp_crick:
        for end in tp_crick[chrom]:

            if end + 2 in fp_watson[chrom]:
                partner = fp_watson[chrom][end+2]
            else:
                partner = 0
            asym_score = (tp_crick[chrom][end] + pc) / (partner + pc)
            print("{}\t{}\t{}".format(chrom, end, asym_score))
    
    
def main():
    desc = '''Report the asymmetry scores for siRNAs. The default is to traverse
    all 3\' ends on the Watson strand and to look for the corresponding 5' ends
    on the Crick strand. You can reverse this behavior (i.e. traverse all 3'
    ends on the Crick strand) by togging --reverse
    '''
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', "--input", type=str, required=True,
                        help='an integer for the accumulator')
    parser.add_argument('--reverse', action='store_true',
                        default=False,
                        help='Use the the 3\' end on the \
                        Crick strand as the reference')

    args = parser.parse_args()

    # fn = "hi5.unox.tncl.alphanodavirus.v0.bed2"
    fn = args.input
    rev = args.reverse

    fh = open(fn)
    # Five prime end and three prime end
    # Mested dict: chrom, position as three levels of keys
    fp_watson = {}
    fp_crick = {}

    tp_watson = {}
    tp_crick = {}

    for line in fh:
        line = line.strip()
        col = line.split()
        chrom = col[0]
        start = int(col[1])
        # Note that this is the exact position of the 3prime end
        end = int(col[2]) - 1
        copy = int(col[3])
        ntm = int(col[4])
        strand = col[5]
        seq = col[6]

        # Make sure all 4 dicts have the same set of chroms
        if chrom not in fp_watson:
            fp_watson[chrom] = {}
            tp_watson[chrom] = {}
            fp_crick[chrom] = {}
            tp_crick[chrom] = {}
            
        if strand == "+":
            if start not in fp_watson[chrom]:
                fp_watson[chrom][start] = 0
            if end not in tp_watson[chrom]:
                tp_watson[chrom][end] = 0
            fp_watson[chrom][start] += copy / ntm
            tp_watson[chrom][end] += copy / ntm
        else:
            if start not in fp_crick[chrom]:
                fp_crick[chrom][start] = 0
            if end not in tp_crick[chrom]:
                tp_crick[chrom][end] = 0
            fp_crick[chrom][start] += copy / ntm
            tp_crick[chrom][end] += copy / ntm

    if not rev:
        sys.stderr.write("Traverse all 3\' ends on the Watson strand\n")
        loading_asym(fp_watson, fp_crick, tp_watson, tp_crick)
    else:
        sys.stderr.write("Traverse all 3\' ends on the Crick strand\n")
        loading_asym_rev(fp_watson, fp_crick, tp_watson, tp_crick)
        

if __name__ == "__main__":
    main()
