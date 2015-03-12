#!/usr/bin/env Rscript

## This script assume an input file that has three columns, with the 1st
## column specifying the coordinate on the gene, 2nd column specifying
## signals on the plus strand and 3rd column specifying the signals on
## the minus strand

function plot_signal_on_gene(a, start, end, tt <- my.title) {
    a <- read.table(input)
    if(ncol(a) == 2) {
        write("The input file only has 2 columns. I will fall back to plotting only the plus strand", stderr())
    } else if(ncol(a)==3) {
        write("The input file has 3 columns.", stderr())
        plot(
        }
    }
}
    
args <- commandArgs(TRUE)
input <- args[1]
start <- args[2]
end <- args[3]
output <- args[4]

bn <- basename(input)
my.title <- paste("Signals on ", bn, sep="")
df <- read.table(input)

