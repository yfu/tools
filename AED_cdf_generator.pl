#!/usr/bin/perl -w 
use strict;
use FileHandle;
use Getopt::Std;
use vars qw($opt_b);
getopts('b:');

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\nAED_cdf_generator.pl: returns values for a cumulative fraction plot based on AED.

USAGE: AED_cdf_generator.pl -b <bin size> <maker1.gff ,maker2.gff ...>\n

Options -b <number between 0 and 1> sets the x axis intervals. 0.025 works well.\n\n";



die($usage) unless $ARGV[0];

#my $FILE = $ARGV[0];
my $BIN_UNIT = $opt_b;
my %HASH;
my $BIN_S = 0;

set_bins();
foreach my $F (@ARGV){
    parse($F);
}
report();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    my @files = sort @ARGV;
    print "AED\t";
    print join("\t", @files);

    foreach my $x (sort keys (%HASH)){
	print "\n$x";
	foreach my $f (sort keys (%{$HASH{$x}})){
	    my $total = $HASH{'1.00'}{$f};
	    my $frac_below = $HASH{$x}{$f}/$total;
	    # more rounding to deal with floating point stuff
	    my $rfrac_below = sprintf("%.3f", $frac_below);
	    print "\t$rfrac_below";
	}
    }
    print "\n";
}
#-----------------------------------------------------------------------------
sub set_bins{

    while ($BIN_S <= 1){
        #I used sprintf to deal with floating point issue in perl
	my $x = sprintf("%.2f", $BIN_S);
	foreach my $file (@ARGV){
	    $HASH{$x}{$file} = 0;
	}
	$BIN_S = $x +$BIN_UNIT;
    }
}
#-----------------------------------------------------------------------------
sub parse {
    my $file = shift;	
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	if ($line =~ /_AED=(\d\.\d+?)\;/){
	    load($1, $file);
	}
	else {
	    next;
	}
    }
    
    $fh->close();
}
#-----------------------------------------------------------------------------
sub load{
    my $value = shift;
    my $file  = shift;

    #go through each fo the values in the bin hash and 
    #add 1 if the value passed into the routine in less 
    #than or equal to the bin value	
    foreach my $x (keys %HASH){
	if ($value <= $x || $value eq $x){
	    $HASH{$x}{$file}++;
	}
    }
}
#-----------------------------------------------------------------------------
