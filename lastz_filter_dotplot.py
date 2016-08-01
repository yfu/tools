#!/usr/bin/env python

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import numpy as np
from itertools import izip


def parse_args():
    parser = argparse.ArgumentParser(description='Filter a dot plot '
                                                 'by LASTZ')
    parser.add_argument('input_file')
    parser.add_argument('min_length', type=int,
                        help='the minimal length of an alignment')
    parser.add_argument('output_file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    # read alignments from the input dot plot file
    with open(args.input_file) as input_file:
        lines = input_file.readlines()
    header = lines[0]
    lines = map(str.split, lines)
    start_positions = np.array(map(lambda x: map(int, x), lines[1::3]))
    end_positions = np.array(map(lambda x: map(int, x), lines[2::3]))

    # for each alignment, determine its length as the average value
    # of the aligned query and target regions
    alignment_len = np.mean(np.abs(start_positions - end_positions),
                            axis=1)

    # filter the alignments
    filtered_alignments = []
    indices = alignment_len > args.min_length
    for i, j in izip(list(start_positions[indices, ]),
                     list(end_positions[indices, ])):
        filtered_alignments.append(map(str, list(i)))
        filtered_alignments.append(map(str, list(j)))
        filtered_alignments.append(('NA', 'NA'))

    # write to the specified output file
    with open(args.output_file, 'w') as output_file:
        output_file.write(header)
        for i in filtered_alignments:
            output_file.write('\t'.join(i) + '\n')
