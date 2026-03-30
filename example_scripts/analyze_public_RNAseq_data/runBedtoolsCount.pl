#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

#### usage - first file is the bed file specifying regions I want to get counts for, the rest of the files are bam files with read mappings
### before Feb 29 2016 this script was calling coverageBed, but now, I use bedtools multicov, e.g. 
# bedtools multicov -split -bams ../*cutadapt.filt.gsnap.sorted.bam -bed flyBaseGene.bed.uniqueID > flyBaseGene.bed.uniqueID.multicov)
## with flyBaseGene.bed.uniqueID I did check that it gave me more or less the same results as before (yes, very similar, but not identical)
## note that now I get just a single output file, rather than one per bam file
## records names of bam files used to run it, and add those as a header to the output file

#my $useSplitOption = 0;  ### for unspliced alignments
my $useSplitOption = 1;  ### for RNA-seq

my $outfileTag = ""; ## this tag will appear in output file names

### set this if on rhino and want to put out one job to each node
my $use_sbatch = 0;
my $walltime = "0-6";

#my $stripTextFromBamFileNameWhenMakingHeader = "_001.cutadapt.filt.gsnap.sorted.bam";
my $stripTextFromBamFileNameWhenMakingHeader = "_STAR.bam";
#my $stripTextFromBamFileNameWhenMakingHeader = ".gsnap.sorted.bam";
#my $stripTextFromBamFileNameWhenMakingHeader = ".first36bp.sorted.bam";
#my $stripTextFromBamFileNameWhenMakingHeader = "";

####### $bedType will be used to construct the beginning of the header (the rest will show which bam files were counted)
my $bedType = "bed12";
my $headerStem = "";
if ($bedType eq "bed12") {
    $headerStem = "chr\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
}
if ($bedType eq "bed12plusID") { # bed12plusID = bed12 with uniqueID added
    $headerStem = "chr\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\tID";
}
if ($bedType eq "bed6") { $headerStem = "chr\tstart\tend\tname\tscore\tstrand"; }
if ($bedType eq "bed3") { $headerStem = "chr\tstart\tend"; }
if ($bedType eq "bed4") { $headerStem = "chr\tstart\tend\tname"; }
if ($bedType eq "gff3") { 
    $headerStem = "chr\tgenome\tfeatureType\tstart\tend\tscore\tstrand\tframe\tinformation";
}
if ($headerStem eq "") {
    die "\n\nTerminating - did not recognize the bedType you supplied. Accepted options are bed12, bed12plusID, bed6, bed3, bed4, gff3\n\n";
}
## get any non-default options from commandline
GetOptions("split=i"      => \$useSplitOption,  ## split reads?
           "outfileTag=s" => \$outfileTag,
           "sbatch=i"     => \$use_sbatch,        ## default 0
           "walltime=s"   => \$walltime,       ## default 0-6 (6hrs)
           "stripText=s"  => \$stripTextFromBamFileNameWhenMakingHeader,
           "bedType=s"    => \$bedType
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";


###########################

if (@ARGV < 2) { die "\n\nterminating - please provide (1) bed file specifying regions I want to get counts for (2) bam file(s) with read mappings to get the counts from\n\n"; }

if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

my $bedfile = shift @ARGV;
my @bamfiles = @ARGV;
if (!-e $bedfile) {die "\n\nterminating bedfile $bedfile does not exist\n\n";}

### get beginning of all outfile names
my $outfileStem = $bedfile;
if ($outfileStem =~ m/\//) { $outfileStem = (split /\//, $outfileStem)[-1]; }
$outfileStem .= ".counts";

if ($useSplitOption == 1) { 
    $outfileStem .= ".split"; 
}
if ($outfileTag ne "") {
    $outfileStem .= ".$outfileTag"; 
}

my $shellScript = "$outfileStem.sh";
my $outfile = "$outfileStem.multicov";
my $outfile2 = "$outfile.temp";

my $options = "";
if ($useSplitOption == 1) { $options .= " -split "; }

### go through the bam files, check each one exists, and construct the header for the output file
foreach my $bamfile (@ARGV) {
    if (!-e $bamfile) {die "\n\nterminating bamfile $bamfile does not exist\n\n";}
    my $shortBamFile = $bamfile;
    if ($shortBamFile =~ m/\//) {
        $shortBamFile = (split /\//, $shortBamFile)[-1];
    }
    $shortBamFile =~ s/$stripTextFromBamFileNameWhenMakingHeader//;
    $headerStem .= "\t$shortBamFile";
}
my $tempFile = "$bedfile.tempheader.txt";
if ($tempFile =~ m/\//) {
    $tempFile = (split /\//, $tempFile)[-1];
}
open (HEADER, "> $tempFile");
print HEADER "$headerStem\n";
close HEADER;

if (-e $outfile) {
    print "    Skipping, as outfile $outfile exists already\n";
    next;
} else {

    open (SH, "> $shellScript");
    print SH "#!/bin/bash\n\n";
    print SH "source /app/lmod/lmod/init/profile\n";
    print SH "module purge\n";
    print SH "module load BEDTools/2.31.0-GCC-12.3.0\n\n";

    ## bedtools multicov
    print SH "bedtools multicov \\\n"; 
    if ($useSplitOption == 1) { print SH "    -split \\\n"; }
    print SH "    -bams";
    foreach my $bam (@ARGV) {
        print SH " $bam";
    }
    print SH " \\\n";
    print SH "    -bed $bedfile > $outfile2\n\n";
    
    ## add header to bedtools output
    print SH "cat $tempFile $outfile2 > $outfile\n\n";
    
    ## tidy up
    print SH "rm $tempFile $outfile2\n\n";
    print SH "module purge\n\n";
    close SH;
    
    ## run the shellScript
    my $run_command;
    if ($use_sbatch == 1) { 
        my $time = "";
        $run_command = "sbatch $time -t $walltime --job-name=bedToolsMulticov $shellScript"; 
    }  else {
        $run_command = "bash $shellScript"; 
    }
    system($run_command);
}

