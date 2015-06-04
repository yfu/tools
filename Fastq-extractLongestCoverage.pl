#!/usr/bin/perl
# From http://wgs-assembler.sourceforge.net/wiki/index.php/Fastq-extractLongestCoverage.pl
use strict;

#  Two pass

die "usage: $0 genomeSize coverage reads.fastq > longest.fastq\n"  if (scalar(@ARGV) != 3);

my $genomeSize = $ARGV[0];
my $coverage   = $ARGV[1];
my $readsIn    = $ARGV[2];

die "Invalid genome size '$genomeSize'\n"  if ($genomeSize == 0);
die "Invalid coverage '$coverage'\n"       if ($coverage == 0);
die "Invalid reads input '$readsIn'\n"     if (! -e $readsIn);

my @lengths;

open(F, "< $readsIn");
while (!eof(F)) {
    my $a = <F>;
    my $b = <F>;
    my $c = <F>;
    my $d = <F>;
    my $l = length($b) - 1;

    push @lengths, $l;
}
close(F);

@lengths = sort { $b <=> $a } @lengths;

my $sumLength = 0;
my $minLength = 0;
my $tgtLength = $genomeSize * $coverage;

my $numReads  = 0;

foreach my $l (@lengths) {
    $sumLength += $l;
    $numReads  += 1;

    if ($sumLength >= $tgtLength) {
        $minLength = $l;
        last;
    }
}

print STDERR "sumLength = $sumLength\n";
print STDERR "minLength = $minLength (with $numReads reads)\n";
print STDERR "tgtLength = $tgtLength\n";

open(F, "< $readsIn");
while (!eof(F)) {
    my $a = <F>;
    my $b = <F>;
    my $c = <F>;
    my $d = <F>;
    my $l = length($b) - 1;

    if ($l >= $minLength) {
        print $a;
        print $b;
        print $c;
        print $d;
    }
}
close(F);
