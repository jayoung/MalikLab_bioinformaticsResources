# Deep sequencing tools

WARNING - MANY OF THESE NOTES ARE VERY OUTDATED!

## general notes

keep an eye on disk space - these will generate large files. Some of the intermediate files can be deleted.

show file size(s) in Mb:
`du -sm filename(s)`

show space left on the disk you're working on
`df -k | grep 'jayoung'`

# quality filtering, trimming, etc

## FASTQC
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

this command without any arguments opens up a GUI interface - I haven't tried it, but it might be useful:
`fastqc`

show usage and arguments:
`fastqc -help`


example, running on a single fastq.gz file:
```
fastqc --outdir fastqc_output/Sample_1_Control_DNA/R2 --format fastq --threads 4 --contaminants /fh/fast/malik_h/user/jayoung/general_notes/NextGenSeqStuff/adapterSequences/variousAdaptersBothStrands.fa.txt ../FASTQ_files/Sample_1_Control_DNA/mySample1_1.fq.gz
```

example, running on all the R2 files in a specified directory to create a single output for each sample (the –casava option tells it to use files names to figure out which fastq files belong to which sample):
```
fastqc --casava --outdir Sample_1_Control_DNA/R2 --format fastq --threads 3 --contaminants /fh/fast/malik_h/user/jayoung/general_notes/NextGenSeqStuff/adapterSequences/variousAdaptersBothStrands.fa.txt ../FASTQ_files/Sample_1_Control_DNA/*R2*gz
```

## filter out read pairs that failed the Illumina chastity filter

Usually the files we get from the sequencing facility ALREADY have failed reads removed. But not always.

I wrote my own R script to do this, as I couldn't find an obvious application I could download to do it. 

The actual work is done by an R script, found here:
~/malik_lab_shared/Rscripts/filterPairedEndReadsIlluminaFilter.R

If you wanted to run that R script on a single pair of reads, you could do this:
```
/fh/fast/malik_h/grp/malik_lab_shared/bin_linux_gizmo/Rdevel/Rscript ~/malik_lab_shared/Rscripts/filterPairedEndReadsIlluminaFilter.R myR1file.fastq.gz myR2file.fastq.gz > & filterPairedEndReadsIlluminaFilter.Rout
```

But if you have multiple pairs of files to filter, it's easiest to call the R script using a perl script, which I will send by email - put this in the bin directory within your home directory. To use it, first open the script in a text editor to make sure the options at the top are set how you'd like them (e.g. whether or not we use sbatch to run in batch mode on the cluster), and then run the script like this:
```
filterPairedEndReadsIlluminaFilterCallRscript.pl my_R1_fastqFile(s)
```

Note - only specify the read1 files as input to the perl script - the script will figure out what the corresponding read2 files are called.  It'll work fine if the files have "_R1" in the name (as well as "fastq" or "fq") - if the file naming convention is different, we'll need minor edits to the perl script so that it can figure out what the reverse reads are called.

Outputs: 
(a) the filtered fastq file names end in something like "filt.fastq.gz", depending on the name of your input file
(b) you'll see some log files appear (names end "IlluminaFilt.Rout") - the end of this file will tell you how many read pairs were kept/filtered.  
(c) you'll also see some slurm log files (e.g. slurm-15090680.out).
If something goes wrong, error messages in files (b) and/or (c) should help figure it out.

## trim off adapters and poor quality sequence, using cutadapt program
https://cutadapt.readthedocs.org/en/latest/guide.html

Show usage and arguments:
cutadapt -help

example commands to trim paired end reads:
- trims based on quality (-q 10) and two adapter sequences that might be on the 3; end
- throws out read pairs less than 20bp long after trimming

first trim based on the forward reads, storing fastq output in temporary files, storing screen output in the cutadaptinfo.txt file:
cutadapt -q 10 a1=adapter1 -a a2=adapter2 --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq reads.1.fastq reads.2.fastq > cutadaptinfo.txt

then trim based on the reverse reads, using the temporary files as input, concatenating screen output to the cutadaptinfo.txt file:
cutadapt -q 15 a1=adapter1 -a a2=adapter2 --minimum-length 20 --paired-output trimmed.1.fastq -o trimmed.2.fastq tmp.2.fastq tmp.1.fastq >> cutadaptinfo.txt

perhaps use my script to call cutadapt more easily on multiple files, e.g.
runCutadapt.bioperl *R1*.fastq.gz

## trim exact basepair amounts off fastq reads using fastx_trimmer
http://hannonlab.cshl.edu/fastx_toolkit/commandline.html

Show usage and arguments:
fastx_trimmer -help

e.g. to trim off the first 10bp of the read (keeps 11-end):
fastx_trimmer -f 11 -i SRR826809_1.fastq -z -o SRR826809_1.trim11toEnd.fastq.gz

or, (I think) using a gzipped input file:
zcat SRR826809_1.fastq.gz | fastx_trimmer -f 11 -z -o SRR826809_1.trim11toEnd.fastq.gz


# Mapping Illumina reads

## BWA

[BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)

First, index the genome:
```
module load BWA/0.7.17-GCC-10.2.0
bwa index refGenome.fa
```

Then map. An example shell script for paired end reads is given [here](../example_scripts/bwa_mem_paired.sh)

# using PacBio reads

Aligners: blasr, minimap2, others

# De novo assemblies

## K-mer analysis using jellyfish (pre-assembly analysis)

example on some tetrahymena data:
1. `jellyfish count`
2. `jellyfish histo`

`-m` specifies kmer size
`--bf-size`:  Bloom filter does something to help reduce memory requirements and output file sizes
`-s` : "The size parameter is an indication of the number k-mers that will be stored in the hash. For sequencing reads, one this size should be the size of the genome plus the k-mers generated by sequencing errors."
`-t` num threads/CPUs
`-C` counts "canonical" kmers (??)
`-o` output file name

```
jellyfish count -m 31 -s 200M --bf-size 10G -t 10 -C -o LGCallFiles.31mercounts.BloomFilter.jf <(zcat ../FASTQ_files_3/LGC_FK01_TPB1hp+p-mixGenomicDNA/LGC_FK01_349bp_CGATGT_L001_R1.fastq.gz)  <(zcat ../FASTQ_files_3/LGC_FK01_TPB1hp+p-mixGenomicDNA/LGC_FK01_349bp_CGATGT_L001_R2.fastq.gz)  <(zcat ../FASTQ_files_3/LGC_FK01_TPB1hp+p-mixGenomicDNA/LGC_FK01_349bp_CGATGT_L002_R1.fastq.gz)  <(zcat ../FASTQ_files_3/LGC_FK01_TPB1hp+p-mixGenomicDNA/LGC_FK01_349bp_CGATGT_L002_R2.fastq.gz) 
     ## used ~ 14 Gb of memory (of 48 Gb total on the machine). Took ~9mins. Output is ~1.5Gb file
jellyfish histo LGCallFiles.31mercounts.BloomFilter.jf > LGCallFiles.31mercounts.BloomFilter.jf.histo.txt
```

use R script to look at tab-delimited results.  Can then calculate coverage and genome size

## De novo assembly using SOAPdenovo2 + GapCloser 

make a folder for your assembly, and work in it, e.g.:
```
mkdir mySOAPassembly1
cd mySOAPassembly1
```
1. set up a config file, calling it whatever you want. Here we copy the example config file to a new file called soap.config and edit as appropriate
```
cp malik_lab_shared/help/SOAP/SOAPdenovo2/example.config ./ soap.config
```
and edit that file - tell it where to find reads, an estimate of library insert size, maximum read length. Lisa and I found for the Bachtrog data that the library insert size parameter made a big difference (500bp, as stated in their methods section, gave a much better assembly than the default of 200bp).  Might even be worth making an initial assembly, mapping reads back to it, calculating average insert size and re-assembling.

2. you'll probably want to run this on a compute node of your own, so you can either start up a grabfullnode session, or you can run this via sbatch:   `sbatch --wrap="myCommand"`
For some datasets, I found SOAP crashed on the regular compute nodes, but worked when I ran it on a "largenode" machine with more memory.

3. run SOAPdenovo2, using 12 threads, kmer size=31, putting screen outputs into files called SOAPdenovo31mers.ass.log and SOAPdenovo31mers.ass.err:
```
(SOAPdenovo-63mer all -s both.config -K 31 -o SOAPdenovo31mers -p 12 > SOAPdenovo31mers.ass.log) >& SOAPdenovo31mers.ass.err
```

4. use GapCloser on the resulting scaffolds file, also using 12 threads
```
GapCloser -b both.config -a SOAPdenovo31mers.scafSeq -o SOAPdenovo31mersGapsClosed -t 12 > SOAPdenovo31mersGapsClosed.stderr
```

5. and then you probably want to format gap-closed scaffolds for blast:
```
makeblastdb -in SOAPdenovo31mersGapsClosed -dbtype nucl -parse_seqids
```

perhaps also check out my script to run all these steps on the compute cluster:
`runSOAP.pl`

SOAP mostly ran OK for me, but sometimes crashed, or ran for a long time without getting anywhere. It does not give helpful error messages. In some cases I fixed that by (a) using a larger memory machine (b) instead of using multiple split fastq files as input, I combined all forward reads into one file, and all reverse reads into a second file.

note for SOAP: "If you're assembling a high heterozygosity genome, you can adjust the -M parameter higher. And if you're assembling a high repeat genome, you can try the -R and -F parameter."
Other assemblers I've read about that might be worth trying, if SOAP doesn't do well: Cortex, Platanus (they claim to do better with heterozygous genomes)

## De novo assembly using platanus (should be SNP-tolerant)

It seems to run fine on a "regular" gizmod node (do not need to wait in the "largenode" queue).  It runs in reasonable time (< 1hr for the three steps)

some documentation in ~/malik_lab_shared/help/platanus/README.txt

There are three steps: assemble, scaffold and gap_close.  Get help on each step:
```
platanus assemble -help
platanus scaffold -help
platanus gap_close -help
```
preparing fastq files: platanus cannot seem to deal with fastq.gz files, so we need to gunzip them.  I also found it simplest to combine all forward reads from one library into a single file (same for reverse reads), e.g.
```
cat *R1*gz > combined_forwardReads.fastq.gz
gunzip combined_forwardReads.fastq.gz

cat *R2*gz > combined_reverseReads.fastq.gz
gunzip combined_reverseReads.fastq.gz
```
make sure the name of the output file for "cat" command does not fit the pattern you're specifying for the input files (e.g. *R1*gz). You can rename the output later if you like.

Running platanus:
Here's how I ran platanus on a single library:
1. assemble
```
platanus assemble -t 12 -o N1C59_plat1 -f /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R1.filt.cutadapt.fastq /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R2.filt.cutadapt.trim7toEnd.fastq 2> ass_log.txt
```
`-t` = threads
`-o` = output file prefix
it allows itself 16Gb of memory for making kmer distribution - this is tunable using the -m parameter

2. scaffold
```
platanus scaffold -t 12 -o N1C59_plat1 -c N1C59_plat1_contig.fa -b N1C59_plat1_contigBubble.fa -IP1 /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R1.filt.cutadapt.fastq /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R2.filt.cutadapt.trim7toEnd.fastq  2> sca_log.txt
```
3. close gaps
```
platanus gap_close -t 12 -o N1C59_plat1 -c N1C59_plat1_scaffold.fa -IP1 /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R1.filt.cutadapt.fastq /fh/fast/malik_h/user/jayoung/michael_ailion_worms/fastqFilesFiltered/N1C59/combinedFiles/N1C59_R2.filt.cutadapt.trim7toEnd.fastq 2>gap_log.txt
```
and the output assembly is called JU1825_plat1_gapClosed.fa

Can specify multiple libraries in steps 2 and 3 as follows:
`-IP1 forward_lib1.fq reverse_lib1.fq  -IP2 forward_lib2.fq reverse_lib2.fq -IP3 forward_lib3.fq reverse_lib3.fq`
In step 1 I think you just list all fastq files in any order (I think it does not use pairing or library information)

## De novo assembly using velvet (and oases for transcriptomes)
https://www.ebi.ac.uk/~zerbino/velvet/

transcriptome assembly pipeline:  1. velveth 2. velvetg 3. oases
genomic seq assembly pipeline:  1. velveth 2. velvetg 

show usage and arguments:
```
velveth
velvetg
oases
```
when Lisa and I assembled the Bachtrog data, velvet did not perform nearly as well as SOAP.

example, using paired ends, kmer length 21bp:
1. velveth, putting output in directory called velvetAssembly
```
velveth velvetAssembly 21 -shortPaired -fastq.gz -separate SRR826657_1.fastq.gz SRR826657_2.fastq.gz
```
2. velvetg, supplying library insert size of 200bp:
```
velvetg velvetAssembly -ins_length 200 -read_trkg yes
```
3. (for transcriptomes only):
```
oases velvetAssembly -ins_length 200
```
4. maybe check number of contigs/transcripts in the output:
```
grep '>' transcripts.fa | wc
grep '>' contigs.fa  | wc
```

# Getting assembly statistics using QUAST 
http://bioinf.spbau.ru/quast

1. make a directory where you will put the results, e.g. 
```
mkdir SOAPdenovo31mersGapsClosed_quastResults
```
2. run the QUAST program: (--threads is number of CPUs to use: default is 8 or 12 so be make sure you ask for appropriate resources)
```
quast.py --threads 3 -o SOAPdenovo31mersGapsClosed_quastResults SOAPdenovo31mersGapsClosed
```
To compare multiple assemblies, you can run quast with several input files, e.g. for two assemblies I had in different subdirectories called BII and LGC_FK01:
```
quast.py --threads 3 -o QUASTtryTwoAssemblies BII/SOAPdenovo31mers.scafSeq.GapsClosed LGC_FK01/SOAPdenovo31mers.scafSeq.GapsClosed
```
Quast will make some recognizable names for your assemblies, but they might end up being really long. An easy cheat to get shorter names in the report is to make aliases (links) with names of your choice:
```
cd BII 
ln -s SOAPdenovo31mers.scafSeq.GapsClosed BII.fa
cd ../LGC_FK01/
ln -s SOAPdenovo31mers.scafSeq.GapsClosed LGC_FK01.fa
quast.py --threads 3 -o QUASTtryTwoAssemblies BII/BII.fa LGC_FK01/LGC_FK01.fa
```

# Using BUSCO to assess assembly quality

BUSCO compares a set of conserved genes against a genome (or transcriptome) to assess completeness, duplication levels, etc

http://busco.ezlab.org

(fhcrc account needs to be set up with correct environmental variables:  
`cp /fh/fast/malik_h/grp/malik_lab_shared/help/loginFileTemplates/.gizmocombinedMalikLabrc_bashVersion ~` )

get help: `BUSCO -help`

example command:
```
BUSCO -in myAssembly.fa -o NameOfOutputDir -l ~/malik_lab_shared/lib/BUSCO/arthropoda -m genome -c 4 -sp fly
```
`-c` = number of threads  
`-l` = the BUSCO library of conserved genes we want to use. Choices:  
~/malik_lab_shared/lib/BUSCO/arthropoda  
~/malik_lab_shared/lib/BUSCO/eukaryota  
~/malik_lab_shared/lib/BUSCO/metazoa  
    (other libraries can be downloaded from http://busco.ezlab.org - vertebrates, fungi, bacteria, plants)

`-sp` = Default: generic. Select species for Augustus to use in making gene predictions. Selecting a closely-­related species usually produces better results. Find out the valid options using this command: "augustus --species=help" 
`-m` = mode. Can be genome/ogs/trans  (ogs=gene set ; trans=transcriptome)

output will go in a folder named "run_NameOfOutputDir".  The most useful file has name starting with "short_summary".   

Probably want to remove many of the other files to save disk space. Unless you are going to use BUSCO output to train a gene prediction program, you probably want to do this: 
```
cd run_NameOfOutputDir
rm -r augustus augustus_proteins gb gffs hmmer_output selected
```

See this folder for a help file, example files, and assessments of various published genomes from the BUSCO website:
`~/malik_lab_shared/help/BUSCO/`

# Using deep-seq data from NCBI's short read archive (SRA)

.sra files are large files, and contain the reads in an archived form.  Quality values might be in several formats - be careful if you use fastq files in downstream applications.

Be careful - these commands will produce large output files.

One way to download files is via the ASPERA program on the mac - download it from here and install it:
http://asperasoft.com/software/transfer-clients/connect-web-browser-plug-in/

Use NCBI's web site to browse and select one or more SRA datasets - when you try to download them it will fire up ASPERA, and you will get an .sra file, e.g. SRR826657.sra

put that dataset on fred
make a link (alias) without the sra extension, as blast doesn't like that:
```
ln -s SRR826657.sra ./SRR826657
```
to directly blast the reads in an sra archive (just blastn, can't do tblastn)
```
blastn_vdb -query testRead.fa -db /home/jayoung/traskdata/public_databases/NCBI/SRA/data/SRR826657/SRR826657 -out testRead.fa.blastnSRR826657
```
To convert sra to fastq format so we can do tblastn, or use the reads in other programs
``` 
fastq-dump --readids --split-spot --split-3 --gzip ../SRR826657.sra
```
(the `--split-spot` option ensures that the two reads for each pair are NOT concatenated together, e.g. without it 2x75bp reads would be concatenated to 150bp reads)
(`--split-3` puts the left and right ends in different files, and any unpaired reads in a third file)
(`--readids` ensures the left and right ends have different seq names)

or, can do it by contacting NCBI website, I think, giving just accession: 
```
fastq-dump --readids --split-spot --split-3 --gzip SRR826657
```    
(I think this method might be quite a bit slower than downloading the file via Aspera and transferring it, and deleting it later, even though it seems easier because there are fewer steps)

# To convert fastq.gz files to fasta format
```
zcat SRR826657_1.fastq.gz | seqret stdin SRR826657_1.fa
zcat SRR826657_2.fastq.gz | seqret stdin SRR826657_2.fa
```
to combine forward and reverse reads into one file, for blasting
```
cat SRR826657_1.fa SRR826657_2.fa > SRR826657bothEnds.fa
```
NOTE:  in these files (for paired-end datasets), each seq has an extension .1 or .2.  This is a good thing if we are trying to blast the entire set (names need to be unique).  But this is a bad thing for some downstream uses (e.g. ARC assembly by reduced complexity) - I think they expect the two reads in the pair to have exactly the same name.  I wrote a script to remove those extensions: stripExtensionFromSRAfastqFiles.pl. It is probably slower than it might be.

tidy up, by deleting the uncombined fasta files: 
```
rm SRR826657_1.fa SRR826657_2.fa
```

format for blast
```
makeblastdb -in SRR826657bothEnds.fa -dbtype nucl -parse_seqids
```
For some datasets, SRA stores a bam file representing reads mapped to a reference genome - in those cases we can get just reads in certain regions from that bamfile - I haven't done it yet, but see instructions here:
http://genometoolbox.blogspot.com/2013/06/download-regional-bam-file-from.html

# Repetitive elements

[xTea](https://www.nature.com/articles/s41467-021-24041-8) package (x-Transposable element analyzer), from Peter Park’s lab.  Nat Commun, 2021


# Aligning reads to a repeat consensus sequence

General outline. Details still need to be fleshed out

1. take a bam file of reads mapped to the "repeat library" genome. Use samtools view to extract reads for the repeat of interest, which should pull out all reads that mapped to any copy of the repeat in the reference genome.

2. use bam2fasta.pl to call seqret to convert reads to fasta format

3. use blastn to align reads to a consensus. Use a script I have that converts blastn output to sam output, which can be viewed in IGV.


# RNA-seq mapping using tophat

These notes are old!   I would probably use STAR now or HISAT2.

1. index the reference genome:
```
bowtie2-build fly_Apr2006.fa fly_Apr2006
```
First argument is the name of the genome seq, in fasta format. Second argument is what we're calling the formatted genome.

2. index the transcriptome annotations (optional). Tophat can use existing annotations to help its mapping. I'm using a gtf file of refGene annotations that I got via UCSC's table browser (refGene.gtf). Command:
```
tophat2 --GTF /home/jayoung/malik_lab/public_databases/UCSC/fly_Apr2006/misc_tracks/refGene.gtf --transcriptome-index ./refGene ./fly_Apr2006
```

3. run tophat to map a single set of paired-end reads.  Command:
```
tophat2 -o hup_tophat_out --min-intron-length 30 --solexa-quals --prefilter-multihits --transcriptome-index /home/jayoung/malik_lab/public_databases/UCSC/fly_Apr2006/indexForTophat/refGene/refGene /home/jayoung/malik_lab/public_databases/UCSC/fly_Apr2006/indexForTophat/fly_Apr2006 ../FASTQ_files/hup_R1_001.fastq.gz ../FASTQ_files/hup_R2_001.fastq.gz --num-threads 3
```
These are the options I chose - might be different for other datasets:
- `--solexa-quals` because these reads are in old solexa format
- `--prefilter-multihits` maps first to genome then transcriptome 
- `--transcriptome-index` specifies the annotation to use
- `--min-intron-length 30` Drosophila have quite small introns (median=68bp). By default tophat will ignore donor/acceptor pairs closer than 70 bases apart. 

After mapping, you might want to do some kind of filtering.  e.g. might filter for reads mapping uniquely.

# samtools for statistics on bam files

Very useful to see how well mapping worked: shows how many reads total, how many mapped, how many properly paired, etc.
```
samtools flagstat bamfile > statsoutputfile.txt
```
count how many mappings there were
```
samtools view –c file.bam
```

# counting reads for each gene

count reads mapping to each annotated transcript.  I'm using a bed file of refGene annotations that I got via UCSC's table browser (refGene.bed), together with coverageBed from the bedtools package. It adds an extra column to the bed file that has the number of reads mapping to each transcript.

command:
```
coverageBed  -split  -counts -abam ../hup_tophat.bam -b /home/jayoung/malik_lab/public_databases/UCSC/fly_Apr2006/misc_tracks/refGene.bed > refGene.bed.counts.split.hup_tophat
```
I then use various R/Bioconductor commands to read in the gene annotations and convert those to RPKMs. For that calculation, need length of each transcript, and total number of reads in each dataset (could use overall total, or number of mapped reads).

