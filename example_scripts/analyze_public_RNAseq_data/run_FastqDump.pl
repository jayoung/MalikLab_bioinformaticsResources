#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

##### usage: on the command line, we specify accessions in one of two ways:
# (a) simply specify accession(s) on command line separated by spaces
# (b) supply the name of a single text file containing one or more accessions, one per line. If lines contain other text the script will split on tab and take the first field

#   runFastqDump.pl accession1 accession2
#   runFastqDump.pl fileWithListOfSRArunIDs.txt

# SRR000001 is a small-ish example that can be used to test (3.3m reads, 454 GS FLX data)

### default settings 
# - sets up a shell script called fastqDump_sbatch.sh that I should then run to ACTUALLY start the jobs going ('sbatch fastqDump_sbatch.sh')
# - shell script will use job arrays and throttling, default running 10 jobs at once
# - will get gzipped output

### examples: 
## get a single accession 
#   runFastqDump.pl SRR000001 ; sbatch fastqDump_sbatch.sh

## don't gzip the resulting fastq files:
#   runFastqDump.pl -gzip=0 SRR000001 ; sbatch fastqDump_sbatch.sh

## run fastqDump on command line without using sbatch
#   runFastqDump.pl -use_sbatch=no SRR000001

## set all jobs going at once using sbatch (without this option we will make a shell script and run the jobs as an array - that's better)
#   runFastqDump.pl -use_sbatch=yes SRR000001

## same but supply a tag to be used as the job-name in sbatch, so we can better track if jobs have finished
#   runFastqDump.pl -use_sbatch=yes -job=myFqDump SRR000001

## show commands that WOULD be run, but don't actually run them:
#   runFastqDump.pl -use_sbatch=yes -debug=1 SRR000001

#####  This script takes each accession, checks whether the .sra file is already downloaded, and if not, downloads the sra file using prefetch (to a cache directory specified below. Then we will use fastq-dump to convert .sra file to .fastq file.
## the .sra files will go here. You will need to MANUALLY delete these files when you're done
my $sraCacheDir = "/fh/fast/malik_h/grp/public_databases/NCBI/SRA/cache/sra";


##### July 9 2020 - now using SRA toolkit v 2.
##### Aug 26 2020 - adding option to create a sbatch script, allowing throttling of jobs. this is now my preferred option

###### $use_sbatch: determines how jobs are distributed over nodes:  
## 'no' means run each job on the command line sequentially
## 'yes' means use sbatch to send jobs to cluster nodes
## 'createSbatchScript' (my preferred option) creates a shell script - I can set up job arrays this way (e.g. I set up a batch of lots of jobs and throttle so that we only run 10 at once).  
#my $use_sbatch = "no";
#my $use_sbatch = "yes";
my $use_sbatch = "createSbatchScript";

# these options only used if $use_sbatch = "createSbatchScript"
my $arrayBatchSizeForSbatchScript = 10;
my $walltime = "0-04:00:00"; # 4 hrs

my $debug = 0;

my $gzip = 1;

my $jobName = "fastqDump";

my $prefetchArgs = "";

## the Brawand multispecies SRA files are stored centrally:
#my $sraCacheDir = "/shared/biodata/ngs/Public/SRA/SRP007412";

## get any non-default options from commandline
GetOptions("use_sbatch=s" => \$use_sbatch,
           "arraySize=i"  => \$arrayBatchSizeForSbatchScript,
           "walltime=s"   => \$walltime,
           "debug=i"      => \$debug, 
           "gzip=i"       => \$gzip,
           "prefetch=s"   => \$prefetchArgs,
           "job=s"        => \$jobName
            ) or die "\n\nterminating - unknown option(s) specified on command line\n\n";

#####################
if ($use_sbatch eq "yes") {print "\n\nUsing srun to parallelize\n\n";}

# parse the accessions requested - might be from the command-line, might be in a file of accessions:
my @accessionsRequested;
foreach my $argument (@ARGV) {
    if (-e $argument) {
        print "\n####### opening file $argument to get accessions\n";
        open (IN, "< $argument");
        while (<IN>) {
            my $line = $_; chomp $line;
            if ($line eq "") {next;}
            if ($line =~ m/\t/) {$line = (split /\t/, $line)[0];}
            if ($line =~ m/\s/) {
                die "\n\nterminating - accession has a space in it - that doesn't look right\n\n";
            }
            push @accessionsRequested, $line;
        }
    } else {
        push @accessionsRequested, $argument;
    }
}

## now test @accessionsRequested to see if I ALREADY have fastq files for any of them in the current dir
my @accessionsToGet;
foreach my $accession (@accessionsRequested) {
    my $test = "$accession"."*";
    my @files = glob "$test";
    #print "accession $accession\n";
    if (@files > 0) {
        print "    files already exist - skipping $accession\n";
        next;
    }
    push @accessionsToGet, $accession;
}
my $numRequested = @accessionsRequested;
my $numToGet = @accessionsToGet;
print "####### getting a total of $numToGet accessions (out of $numRequested requested)\n";
if ($numToGet == 0) {
    die "\n\nterminating - looks like there are no remaining accessions to get\n\n";
}


my $gzipOption = "";
if ($gzip == 1) { $gzipOption .= "--gzip "; }

##### set up a sbatch script, perhaps
if ($use_sbatch eq "createSbatchScript") {
    ## it wouldn't make sense for batch size to be bigger than the number of jobs (although it might not matter?)
    if ($arrayBatchSizeForSbatchScript > $numToGet) {
        $arrayBatchSizeForSbatchScript = $numToGet;
    } 
    
    open (OUT, "> fastqDump_sbatch.sh");
    print OUT "#!/bin/bash\n";
    print OUT "#SBATCH --nodes=1\n";
    print OUT "#SBATCH --cpus-per-task=1\n";
    print OUT "#SBATCH --time=$walltime\n"; 
    my $arraySize = $numToGet - 1;
    print OUT "#SBATCH --array=[0-$arraySize]\%$arrayBatchSizeForSbatchScript # how many jobs and batch size\n";
    print OUT "#SBATCH --output=slurm.%J.out\n"; # stderr goes here too

    print OUT "#SBATCH --job-name=\"$jobName\"\n\n";
    
    print OUT "ACCESSIONS=(@accessionsToGet)\n\n";
    print OUT "CACHE=\"$sraCacheDir\"\n";
    print OUT "SINGLE_ACCESSION=\"\${ACCESSIONS[\$SLURM_ARRAY_TASK_ID]}\"\n";
    print OUT "SRAFILE=\"\${CACHE}/\${ACCESSIONS[\$SLURM_ARRAY_TASK_ID]}.sra\"\n\n";

    print OUT "source /app/lmod/lmod/init/profile\n";
    print OUT "module load SRA-Toolkit/3.1.1-gompi-2023b\n\n";

    print OUT "echo \"My accession: \${SINGLE_ACCESSION}\"\n";
    
    ## run prefetch if the sra file does not already exist
    print OUT "if test -f \"\${SRAFILE}\"; then\n";
    print OUT "    echo \"\$SRAFILE already exists\"\n";
    print OUT "else\n";
    print OUT "    echo \"\"\n";
    print OUT "    echo \"running prefetch\"\n";
    print OUT "    prefetch $prefetchArgs \${SINGLE_ACCESSION}\n";
    print OUT "fi\n\n";
    
    ## run fastq-dump on the SRA file
    print OUT "echo \"\"\n";
    print OUT "echo \"running fastq-dump\"\n";
    print OUT "fastq-dump --split-spot --split-e $gzipOption \${SRAFILE}\n";
    close OUT;
    print "set up sbatch script - you should run it now, e.g.\n    sbatch fastqDump_sbatch.sh\n\n";
} else {
    foreach my $accession (@accessionsToGet) {
        print "\n###### working on accession $accession\n";
        
        ## for single end data the outfile name will look like this: SRR1640128.fastq.gz
        ## for paired end data the outfile name will look like this: SRR5363157_1.fastq.gz and SRR5363157_2.fastq.gz  
        
        ## $outfile is what it will be called if it's paired end data. $outfileSingle is what it will be called if it's single-end data. 
        my $outfile = $accession; 
        my $outfileSingle = $accession; 
        $outfile .= "_1";
        $outfile .= ".fastq"; $outfileSingle .= ".fastq";
        if ($gzip == 1) { $outfile .=".gz"; }
        if ((-e $outfile) || (-e $outfileSingle)) {
            if (-e $outfile) {
                print "    skipping $accession - outfile $outfile exists already\n\n";
            }
            if (-e $outfileSingle) {
                print "    skipping $accession - outfile $outfileSingle exists already\n\n";
            }
            next;
        }
        
        ## if needed, get the .sra file
        my $sraFile = "$sraCacheDir/$accession.sra";
        my $command = "";
        if (!-e $sraFile) {
            print "    getting sra file using prefetch\n";
            if ($use_sbatch eq "yes") {
                $command .= "prefetch $accession ; ";
            } else {
                $command .= "prefetch $accession ; ";
            }
        } else {
            print "    already have sra file - skipping prefetch, and going direct to fastq-dump\n";
        }
        
        ## use the .sra file to get the fastq files
        $command .= "fastq-dump --split-spot --split-e $gzipOption ";
        $command .= $sraFile;
        if ($use_sbatch eq "yes") {
            $command = "sbatch --job-name=$jobName --wrap=\"$command\"";
        }
        print "command:\n$command\n\n";
        if ($debug == 0) { system($command); }
    }
    
    if ($use_sbatch eq "yes") {
        print "\n\nSet all jobs going - use sq command to monitor whether there are still any fastq-dump commands running\n\nKeep an eye on memory using df -k | grep 'jayoung'\n\n";
    }
}

