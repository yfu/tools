#!/usr/bin/env python

import fileinput

cur_name = ""                   # The name of the current sequence
cur_seq = ""
for line in fileinput.input():
    line = line.strip()
    if(line[0] == '>'):
        if(cur_name != ""):
            # This is not the first sequence
            print cur_name + "\t" + str(len(cur_seq))
            cur_seq = ""
        cur_name = line[1:].strip()
    else:
        cur_seq += line

# Print out the last sequence length
print cur_name + "\t" + str(len(cur_seq))
