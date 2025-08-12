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
GetOptions("split=i" => \$useSplitOption,  ## split reads?
           "outfileTag=s" => \$outfileTag,
           "sbatch=i" => \$use_sbatch,        ## default 0
           "walltime=s" => \$walltime,       ## default 0-6 (6hrs)
           "stripText=s" => \$stripTextFromBamFileNameWhenMakingHeader,
           "bedType=s" => \$bedType
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";


###########################

if (@ARGV < 2) { die "\n\nterminating - please provide (1) bed file specifying regions I want to get counts for (2) bam file(s) with read mappings to get the counts from\n\n"; }

if ($use_sbatch == 1) {print "\n\nUsing sbatch to parallelize\n\n";}

my $bedfile = shift @ARGV;
my @bamfiles = @ARGV;
if (!-e $bedfile) {die "\n\nterminating bedfile $bedfile does not exist\n\n";}
my $logfile = "$bedfile.counts.Logfile.txt";
my $outfile = "$bedfile.counts.multicov";
if ($outfile =~ m/\//) { $outfile = (split /\//, $outfile)[-1]; }
if ($logfile =~ m/\//) { $logfile = (split /\//, $logfile)[-1]; }

if ($useSplitOption == 1) { 
    $outfile =~ s/counts/counts.split/; 
    $logfile =~ s/counts/counts.split/;
}
if ($outfileTag ne "") {
    $outfile =~ s/counts/counts.$outfileTag/;
    $logfile =~ s/counts/counts.$outfileTag/;
}
my $outfile2 = "$outfile.temp";

my $options = "";
if ($useSplitOption == 1) { $options .= " -split "; }

open (LOG, ">> $logfile");
print LOG "#####################\n";
print LOG scalar localtime, "\n";
print LOG "#####################\n";

### go through the bam files, check each one exists, and construct the header
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

print "  counting reads for regions in $bedfile in these bam files: @ARGV\n\n";
print LOG "  counting reads for regions in $bedfile in these bam files: @ARGV\n\n";
    
if (-e $outfile) {
    print "    Skipping, as outfile $outfile exists already\n";
    print LOG "    Skipping, as outfile $outfile exists already\n";
    next;
} else {
    my $command = "bedtools multicov ";
    if ($useSplitOption == 1) { $command .= "-split "; }
    $command .= "-bams @ARGV -bed $bedfile > $outfile2";
    $command .= " ; cat $tempFile $outfile2 > $outfile";
    $command .= " ; rm $tempFile $outfile2";

    if ($use_sbatch == 1) { 
        my $time = "";
        if ($walltime ne "default") { $time = "-t $walltime"; }
        $command = "sbatch $time --job-name=bedToolsMulticov --wrap=\"$command\""; 
    } 
    print "    Running command:\n    $command\n";
    print LOG "    Running command:\n    $command\n";
    system($command);
}

close LOG;


if ($use_sbatch == 1) {
    print "\n\nSet all jobs going - use sq command to monitor whether there are still any bedtools commands running\n\n";
} else {
    print "\n\nFinished - see logfile $logfile\n\n";
}

