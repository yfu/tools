#!/usr/bin/python

# Compare the SNPs between different lines
# Wirte to STDOUT the output matrix and write to STDERR about the progress
import sys

fh = open(sys.argv[1])
line = fh.readline().strip()
e = line.split(',')
# print e
# ln is drosophila line
chrom, ln, coverage = e[0], e[1:-2], e[-2] # The last one is empty
order = ln # So I can keep the order
snp = {}
for l in ln:
    snp[l] = {}

counter = 0
for line in fh.readlines():
    counter += 1
    if counter % 1000 == 0:
        sys.stderr.write("I have read %d lines...\n" % (counter))

    e = line.split(',')
    loc, base, cov = e[0], e[1:-2], e[-2]
    for i in range(len(order)):
        snp[order[i]][loc] = base[i]

# Pairwise comparison
diff_mat = {}
for l in ln:
    diff_mat[l] = {}
for i in order:
    for j in order:
        # TODO: Save some computation: calculate half of the matrix
        l1 = snp[i]
        l2 = snp[j]
        locs = l1.keys()
        diff = 0
        for l in locs:
            if l1[l] != l2[l]:
                diff += 1

        # print "Finish comparing %s, %s" % (i, j)
        # print "Current diff is %d" % (diff)
        diff_mat[i][j] = diff
        sys.stderr.write("Finished comparing %s and %s.\n" % (i, j))
# print len(diff_mat)
# print len(diff_mat[i])
# Report the matrix
# First row
print " ",
for o in order:
    print "\t", o,
print 
# The rest of the rows
for i in order:
    print i,
    for j in order:
        print "\t", diff_mat[i][j],
    print 
