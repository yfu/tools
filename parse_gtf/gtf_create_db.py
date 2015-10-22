#!/usr/bin/env python
# Given a GTF file, output a GTF database

import gffutils
import sys

# gtf="/data/fuy2/shared/hg19/gencode_r19/gencode.v19.annotation.gtf"
gtf = sys.argv[1]
output = sys.argv[2]

db = gffutils.create_db(gtf, output, disable_infer_transcripts=True, disable_infer_genes=True)
