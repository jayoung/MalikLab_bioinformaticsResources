# Repetitive regions, and how to analyze them 

Gene families, transposons, and other repeat types (e.g. satellites) are difficult to analyze, due to multi-mapping short reads.

## Technologies

Long reads (e.g. PacBio, Nanopore) are more uniquely mapable, but we usually get fewer reads. Short reads (Illumina) are more common for counting-based assays (e.g. RNA-seq with triplicates).

## Mapping tools

### RNA-seq mappers

Ching-Ho thinks STAR might do better with multi-mapping reads

#### STAR


#### HISAT2

### DNA-seq mappers

#### bwa

## Read counting tools

### htseq-count

### bedtools 

bedtools multicov

bedtools coverageBed (I used to use this)