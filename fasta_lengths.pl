#!/usr/bin/env perl

while(<>) {
    chomp;
    if(m/\s/) {
	print("No whitespaces allowed in the fasta headers. Consider using cut -f1 to modify the headers. Bye!\n");
	exit 1;
    }
    if(m/^>(.+)/){$id=$1} else{ $l{$id} += length($_); } 
}

for my $i (keys %l) {
    print $i . "\t" . $l{$i} ."\n" 
}
