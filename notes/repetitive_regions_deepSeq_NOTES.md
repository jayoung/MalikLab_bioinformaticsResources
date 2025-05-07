# Repetitive regions, and how to analyze them 

Gene families, transposons, and other repeat types (e.g. satellites) are difficult to analyze, due to multi-mapping short reads.

## Technologies

Long reads (e.g. PacBio, Nanopore) are more uniquely mapable, but we usually get fewer reads. Short reads (Illumina) are more common for counting-based assays (e.g. RNA-seq with triplicates).

## Mapping tools

Reference genomes - be careful whether you have unmasked, soft-masked (lower case) or hard-masked (NNN) version of the reference genome. Some mappers pay attention to soft-masking (maybe less likely to map there?).

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

# methods used by public expression catalogs

TCGA, Human Protein Atlas, GTEX - Sophie says there is documentation of how they mapped reads and counted.

# tools to check out!


# Helena Reyes Gopar, April 2025

We had a visit in April 2025 from Helena Reyes Gopar hreyesgopar@northwell.edu who is a postdoc who recently started in the lab that wrote Telescope (?) and has developed a similar tool to work on single cell data called Stellarscope.


[Telescope](https://github.com/mlbendall/telescope?tab=readme-ov-file) is for transposon expression analysis (bulk data)

They seem very motivated to get expression estimates of individual repeat instances. Uses some sort of EM to assign read counts to the most likeluy unique copy. Bayesian model-based reassignment.

UMI-based deduplication will be messed up for multimappers because it uses map location as well as the barcode sequence, so I think she talked about ways around that.



# chat with Sophia Kogut, Feb 11, 2025 (Daniel Blanco-Melo lab)

Sophie is mostly interested in HERV expression, and trying to quantify expression of individual insertions (not aggregating by class)

She uses hg19, and custom HERV annotations. The annotation you choose makes a difference, and the publically available ones were not good. 

Sometimes they look at all repeats, not just HERVs, and here they ignore satellite seqs.

They sometimes use STAR for mapping and featureCounts for counting, sometimes Salmon.

STAR - they DO allow multi-mapping, and then in featureCounts (standalone, although it is in Rsubreads package) - there is a fractional counts option in featureCounts.  They normalize counts by repeat length, and (when looking at classes of repeat) by number of repeats in the group.

They sometimes use Salmon for mapping (and it does its own quantification). Sean from Tapscott lab said Salmon and Kallisto give similar results.  Salmon does some sort of graph-based alignment and is less rigid. It does "pseudo-mapping".

For Sita Kugel's data, they looked at all TEs and used STAR. For HERVs alone, they use Salmon.

Hutch bioinformatics core also runs some sort of TE analysis pipeline involving TEtranscripts tool (Molly Hammell), but they don't really care much about the details.

I think she said Telescope is worth checking out, as is TEtranscripts (Molly Hammell)

Possibly useful plot: 
- x-axis is log-fold change
- y-axis is different categories of TEs
- point size reflects p-value of log-fold change

kmers - Erin (Barnett) from their lab has tried kmer counting a bit.