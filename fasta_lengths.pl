#!/usr/bin/env perl

while(<>) {
    chomp; if(m/^>(.+)/){$id=$1} else{ $l{$id} += length($_); } 
}

for my $i (keys %l) {
    print $i . "\t" . $l{$i} ."\n" 
}
