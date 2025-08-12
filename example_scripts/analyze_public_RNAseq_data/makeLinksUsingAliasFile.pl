#!/usr/bin/perl

use warnings;
use strict;

my $accsFile = "../human_spermatogenesis_SRAaccs.txt";
open (IN, "< $accsFile");
while (<IN>) {
    my $line = $_; chomp $line;
    my $acc = (split /\t/, $line)[0];
    my $testFile = "../../human_testis_sperm/fastqFiles/single/$acc.fastq.gz";
    if (-e $testFile) { print "\nlinking $testFile\n";  system("ln -s $testFile ."); }
}
close IN;
