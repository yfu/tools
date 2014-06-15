#!/bin/bash
F=${1%.raw}
adapter=TGGAATTCTCGGGGG
Extract_insert_6mer.pl $1 $adapter > $1.insertsout
inserts2uniqreads.pl $1.insertsout 18 30 > $F.inserts
rm $1.insertsout

