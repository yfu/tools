#!/usr/bin/awk -f

# Some programs generate a fastq file that contains both read1 and read2. This simple AWK script
# can split such a file into two separate files (stdout for r1 and stderr for r2)
# Usage: awk -f separate_r1_r2.awk  Hi5_NoIndex_L006.total.mp.fastq > Hi5_NoIndex_L006.total.mp.r1.fastq 2>Hi5_NoIndex_L006.total.mp.r2.fastq
#
# Author: Yu Fu

BEGIN{r=0}
{
    if(NR%4==1)
    {
	if($2=="1:N:0:") {r=1}
	else if($2=="2:N:0:") {r=2}
	else {print "Error!"}
    }    
    if (r==1) {print $0 > "/dev/stdout"}
    else {print $0 > "/dev/stderr"}
}
    
    
