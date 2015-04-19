#!/usr/bin/perl

use strict;
use warnings;

my $format = "";

# set regular expressions
my $sanger_regexp = qr/[!"#$%&'()*+,-.\/0123456789:]/;
my $solexa_regexp = qr/[\;<=>\?]/;
my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\[\]\^\_\`abcdefgh]/;
my $all_regexp = qr/[\@ABCDEFGHI]/;

# set counters
my $sanger_counter = 0;
my $solexa_counter = 0;
my $solill_counter = 0;


my $i;
while(<>){
    $i++;

    # retrieve qualities
    next unless $i % 4 eq 0;

    #print;
    chomp;

    # check qualities
    if( m/$sanger_regexp/ ){
	$sanger_counter = 1;
	last;
    }
    if( m/$solexa_regexp/ ){
	$solexa_counter = 1;
    }
    if( m/$solill_regexp/ ){
	$solill_counter = 1;
    }
}

# determine format
if( $sanger_counter ){
    $format = "sanger";
}
elsif( !$sanger_counter && $solexa_counter ){
    $format = "solexa";
}
elsif( !$sanger_counter && !$solexa_counter && $solill_counter ){
    $format = "illumina";
}

print "$format\n";
