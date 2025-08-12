#!/usr/bin/perl

use warnings;
use strict;
use File::Spec;

##### count reads in each fastq file specified on command line
##### uses wc on command line. Can cope with gzipped files, using zcat if needed
##### usage    countFastqReads.pl *gz *gzip *bz2 *fastq *fq

open (OUT, "> readcounts.txt");
print OUT "File\tDir\tNumReads\n";
my $total = 0;
foreach my $file (@ARGV) {
    if (!-e $file) {
        die "\n\nterminating - cannot open file $file\n\n";
    }
    print "working on file $file\n";
    
    my ($volume,$dir,$shortfile) = File::Spec->splitpath( $file );
    $dir =~ s/\/$//;
    my $numlines = "NA";
    if ($file =~ m/gz$|gzip$/) {
        $numlines = `zcat $file | wc -l`;
    } 
    if ($file =~ m/bz2/) {
        $numlines = `bzcat $file | wc -l`;
    }
    if ($numlines eq "NA") {
        $numlines = `wc -l $file`;
    }
    $numlines =~ s/^\s+//g;
    $numlines = (split /\s+/, $numlines)[0];
    my $numseqs = $numlines / 4;
    print OUT "$shortfile\t$dir\t$numseqs\n";    
    $total += $numseqs;
}
print OUT "\nTotal\t$total\n";
close OUT;
