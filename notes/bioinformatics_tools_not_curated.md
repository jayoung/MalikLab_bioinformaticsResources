# Search tools

NCBI's Edirect toolkit to help search and extract data from their databases. There's an online [book](https://www.ncbi.nlm.nih.gov/books/NBK179288/) and [cookbook](https://ncbi-hackathons.github.io/EDirectCookbook/)

NCBI's SRAtoolkit

LOGAN search to search SRA database using kmers

LexicMap - if you want to search millions of prokaryotic genomes at once, blast can be too slow. This is a quicker alternative.


# Multiple sequence alignments

## Making multiple sequence alignments

Fast Statistical Alignment
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000392
http://orangutan.math.berkeley.edu/fsa/, http://fsa.sourceforge.net/.

## Alignment masking

There are block-filter methods (remove entire columns from alignment) and segment-filter methods (remove part of whole of individual sequences).  I think also some hybrid methods?

MACSE has some sort of mode where it trims non-homologous fragments

Figure 1 of [the OrthoMam paper](https://academic.oup.com/nar/article/52/D1/D529/7318103) shows the pipeline they used for alignments - they do a lot of filtering. 15,868 alignments from 190 mammalian genomes. Website isn't working well as of Oct 11 2024.

### Filter out entire sequences
- [PhylteR](https://academic.oup.com/mbe/article/40/11/msad234/7330000)

### Block methods:

- trimal - can trim gappy positions out of an alignment (nucleotide or protein). Example usage:
```
trimal -in myAln.fa -out myAln.trimal.fa -automated1 -htmlout myAln.trimal.html -colnumbering > myAln.trimal.retainedColumns.txt
```
Some papers say trimal isn't as good as GUIDANCE and/or ZORRO

Also:
- BMGE
- GBLOCKs - some papers say this isn't as good as GUIDANCE

### Segment methods
- [HMMCleaner](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-019-1350-2)

### Other methods - maybe segment, maybe hybrid, or I just haven't looked

I think 
- `SWAMP`
- GUIDANCE (I used this in the kinetochore pipeline)
- ZORRO (Wu, Chatterji, Eisen, PLoS ONE, 2012)
- PREQUAL
- MACSEv2 "trimNonHomologousFragments" module

## Public alignment/orthology datasets

Ensembl and the BioMart interface

OrthoDB (there is an R package for OrthoDB). Provides estimate of evolutionary rates

[OrthoMaM](https://academic.oup.com/nar/article/52/D1/D529/7318103) - 15,868 alignments from 190 mammalian genomes, should be 1-to-1 othologs. Nice alignment filtering: should be high quality? Should be in-frame.

Zoonomia project: >240 mammals analyzed. Precomputed genome-wide mammalian conservation scores, and codon alignments are here https://zoonomiaproject.org/the-data/ 

[CESAR tool](https://github.com/hillerlab/CESAR2.0) somehow extracts gene/exon alignments from whole genome multiple alignments. 

12-species alignments (10 mammals+chicken+tortoise) (nucleotide and amino acid) for 7445 genes from Shuler & Hagai paper 2022 in teh associated [github repo](https://github.com/galshuler/Evolution_of_host-virus_protein_protein_interactions). I downloaded and unpacked in `~/FH_fast_storage/paml_screen/Adrian_Pelin_vaccinia/references/2022-Shuler/git_repo/Evolution_of_host-virus_protein_protein_interactions`. See notes in `~/FH_fast_storage/paml_screen/Adrian_Pelin_vaccinia/Adrian_Pelin_vaccinia_NOTES.md`

120-mammal alignment described in this 2019 paper by [Hecker and Hiller](https://pubmed.ncbi.nlm.nih.gov/31899510/) and available [in a genome browser](https://genome-public.pks.mpg.de/) and [for download](https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/).

144-vertebrate whole genome alignment described in this 2017 paper by [Sharma and Hiller](https://academic.oup.com/nar/article/45/14/8369/3875570) - data are available [here](https://bds.mpi-cbg.de/hillerlab/144VertebrateAlignment_CESAR/). The CESAR tool takes whole genome alignments and human exon coordinates, and uses that to extract orthologous exons.

TreeFam

[Genomicus browser v93](https://www.genomicus.bio.ens.psl.eu/genomicus-93.01/cgi-bin/search.pl) - uses orthology data from ensembl. Older versions like v93 represent paralogs better. 

COG: NCBI's database of clusters of orthologous genes.  [COG website](https://www.ncbi.nlm.nih.gov/research/cog) allows you to send selected protein sequences from COG to COBALT alignment generator, to tweak alignment parameters if needed, and then to download an alignment.

PhyloPat
Seems to be out-of-date and unmaintained
June 2014: http://www.cmbi.ru.nl/cdd/phylopat/52/
Based on Ensembl 52 (quite old:  Ensembl now up to v75)

[xProtCAS](https://slim.icr.ac.uk/projects/xprotcas) server has precalculated scores for many protein structures (or perhaps just structural modules within many protein structures). Can also view/download the underlying alignments

[Consurf-db](https://consurfdb.tau.ac.il/overview.php) (that data also used/displayed in Proteopedia database)
- per-residue rates for proteins that have structure in PDB. 
- Uses up to 300 homologues, found in UniRef90, with some seq identity/coverage filters 
- evolutionary rates for each site calculated based on tree+alignment using Rate4Site algorithm
- Proteopedia warns that the homolog gathering is automated, and may include proteins of similar ancestry but divergent function. To truly understand which residues are important for a particular function, you'd want to filter input protein alignment based on functional knowledge

[orthogene](https://bioconductor.org/packages/release/bioc/vignettes/orthogene/inst/doc/orthogene.html) R package

eggNOG [database](http://eggnog6.embl.de) and [publication](https://pubmed.ncbi.nlm.nih.gov/36399505/)

phyloP tracks from UCSC genome browser

[OMA orthology database](https://pubmed.ncbi.nlm.nih.gov/29106550/) 

INPARANOID (?) maybe old


# Phylogenies

## Making trees

IQ-Tree - I see this used in several large-scale studies

[Treebest](http://treesoft.sourceforge.net/treebest.shtml) - used in Ensembl pipeline. Combines multiple tree types (nucleotide, protein, dN, dS, species tree), tries to reduce number of mismatches to species tree, infers duplication, deletion, etc.

SATe - simultaneous alignment and tree estimation - Liu, Raghavan, Nelesen, Linder and Warnow. Science, 2009


PLATO for detecting variation in evolutionary rates along alignment

## Analyzing trees

NOTUNG for tree reconciliation

CAFE for tree reconciliation (David Liberles says may be too simplistic)

David Liberles: expectation maximization for missing data

OrthoMCL for clustering genes/orthologs

## Getting species trees

[TimeTree](http://timetree.org)

[phyloT](https://phylot.biobyte.de) for getting accepted species trees

NCBI taxonomy

Open Tree of Life
https://tree.opentreeoflife.org/opentree/argus/opentree9.1@ott93302

Lifemap displays of NCBI taxonomy, or the open tree of life:
http://lifemap.univ-lyon1.fr/

# Molecular evolution

## PhyloSim

PhyloSim (http://bit.ly/o09JXW) is an extensible and feature-rich framework for the Monte Carlo simulation of sequence evolution. Can use many different models (nucleotide, codon, amino acid, indels).

PhyloSim is implemented as an R package and it is available from CRAN (http://bit.ly/n54VGg) and the GitHub repository (http://bit.ly/nFo9rf) which also contains many example scripts implementing simulations of varying complexity.


# Phylogenetic tree display/plotting

interactive Tree Of Life [iTOL](https://itol.embl.de) has some nice-looking plotting tools

[ggtree](https://bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/ggtree.html) for R

# Transposons

TEtools - [paper](https://www.ncbi.nlm.nih.gov/pubmed/27924026) and [repo](https://github.com/l-modolo/TEtools)

TEanalysis

RepEnrich2

Molly Hammell (CSHL) may have some tools available, but I see nothing online as of June 2014

Ting Wang may have some tools available, but nothing obvious from his webpage as of June 2014. I think there are some publications I whould look at instead. http://wang.wustl.edu/software

Tedna: a transposable element de novo assembler http://www.ncbi.nlm.nih.gov/pubmed/24894500?dopt=Abstract

RepeatModeler/RECON: de-novo repeat identification, classification and library building software suite. www.repeatmasker.org/RepeatModeler.html

Annotating repeats in a novel species:  EDTA tool, 2019
https://www.ncbi.nlm.nih.gov/pubmed/31843001?dopt=Abstract

TypeTE - Genotyping known polymorphic TE insertion loci:  https://www.ncbi.nlm.nih.gov/pubmed/32067044?dopt=Abstract

xTea (x-Transposable element analyzer) Peter Park lab 2021 https://pubmed.ncbi.nlm.nih.gov/34158502/

TEHub resources Travis Wheeler https://pubmed.ncbi.nlm.nih.gov/34154643/


# User-friendly online toolsets

Galaxy
- Main public instance
- Hutch instance
- various other instances

EDGE - seems to be mostly for bacterial / metagenomics analyses?

# Deep sequencing

GC normalization: John Huddleston advised using mrFAST or mrsFAST for read mapping and then mrcanavar for GC normalization and copy number estimation. But that was a while ago (<2017?)

## RNA-seq
Ching-Ho likes HISAT2 (new version of tophat) for RNA-seq, partly because it has fewer options than STAR, but he also thinks STAR might do better with multi-mapping genes.  HISAT2 may be able to directly take SRA accessions and map without downloading reads.

Ching-Ho also notes that some RNA-seq datasets have high rRNA contamination. This could be affecting my RPKMs quite a lot.

# Datasets/resources

## Genome browsers

UCSC
for custom tracks, hubs, etc, we can use a directory on fast that is accessible from beyond the Hutch firewall: /fh/fast/malik_h/pub/
see the notes file in that directory for more detail on how to host a hub (malik_pubDir_NOTES_JY.txt)

Ensembl

## Species-specific databases

Flybase, Wormbase, SGD

Virus Pathogen Resource (ViPR) is now part of [BV-BRC](https://www.bv-brc.org), the bacterial and viral bioinformatics resource center




## Methods for conservation / evolutionary rates

### Nucleotide

dN/dS type methods

genome-based methods (phyloP etc)

[SLR](https://www.ebi.ac.uk/research/goldman/software/slr/)

Selecton, 2007 - [publication](https://pubmed.ncbi.nlm.nih.gov/17586822/) and [website](http://selecton.bioinfo.tau.ac.il). Website not responding on 10/11/2024. I don't see a way to download the algorithm. Looking at the paper, it looks like it uses various evolutionary models (e.g. codeml M7/M8 etc, or MEC, which takes protein models like JTT and converts them to codon models). Does some kind of hypothesis testing, chooses the best-fitting model.  Also uses a Bayesian method to estimate sitewise dN/dS - uses 8 discrete rate categories and does something with posteriors to get the estimate.

### Protein

[Evorator](https://evorator.tau.ac.il) - website not responding on 10/4/24

[Rate4Site](https://www.tau.ac.il/~itaymay/cp/rate4site.html) (ref [Pupko et al, 2002](https://pubmed.ncbi.nlm.nih.gov/12169533/)) - given a multiple sequence alignment, it can generate a tree (or you can provide one) and it calculates the relative evolutionary rate at each site.

## Transcriptomics

Multispecies developmental timecourse [browser](https://apps.kaessmannlab.org/evodevoapp/).
- human, rhesus macaque, mouse, rat, rabbit, opossum and chicken
- seven organs (cerebrum, cerebellum, heart, kidney, liver, ovary and testis)
- developmental time points from early organogenesis to adulthood

## Ribosomal profiling (ribo-seq)

https://github.com/nzhang89/RiboSeeker

This paper PROBABLY has some sort of accessible database associated [High-quality peptide evidence for annotating non-canonical open reading frames as human proteins](https://www.biorxiv.org/content/10.1101/2024.09.09.612016v1)

## Single cell data

Brotman Baty's [Descartes](https://descartes.brotmanbaty.org) - gateway to single cell atlas for several species

[This paper](https://www.nature.com/articles/s41586-024-08411-y) from Aviv Regev's lab shodescribes a tool called `SCimilarity` that allows you to take a single cell profile query and search all known existing single cell data to find datasets that contain cells with matching profiles.

[This 2024 perspective paper](https://www.nature.com/articles/s41586-024-08338-4) from Aviv Regev and colleagues might be a good starting point for learning about uses of single cell data


## Variation 

Human: 1000 genomes [selection browser](http://hsb.upf.edu/)

Drosophila: [popfly](http://popfly.uab.cat/) (has data to allow M-K tests)


# Variation

SHAPEIT for phasing haplotypes

## Databases of variation

### SNPs

[gnomAD](https://gnomad.broadinstitute.org)

[European Variant Archive](https://www.ebi.ac.uk/eva/) - I don't think it shows allele frequencies

Ensembl browser also has table view that allows filtering. [Example](https://useast.ensembl.org/Homo_sapiens/Gene/Variation_Gene/Table?g=ENSG00000171497;r=4:158709127-158723396)  (but in this example, it's missing a variant that's common according to gnomAD)

dbSNP

dbVar for copy number/structural variants

ClinVar for clinically significant variants

UCSC browser

[Human Genetic Variation Database](https://www.hgvd.genome.med.kyoto-u.ac.jp) (Kyoto)

1000 genomes https://www.internationalgenome.org

### Structural variants

[Database of Genomic Variants](https://www.ebi.ac.uk/dgva/) - being phased out, transitioning to [European Variant Archive](https://www.ebi.ac.uk/eva/)

## Variant callers for long reads

Poster at Genome Informatics meeting said Clair3 works really well for nanopore bacterial sequencing


# Protein sequence analysis

for structure, Antoine used UCSF Chimera (Pettersen, et al. 2004).

AFDB - alphaFold protein structure database

BFVD - database of predicted virus structures
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkae1119/7906834#google_vignette

STRINGdb - database of known and predicted protein-protein interactions

[DHR](https://pubmed.ncbi.nlm.nih.gov/39123049/) - dense homolog retrieval. Use cases: 
- distand homology
- when speed is important
Alignment-free. Mysterious encoding of protein seqs, allows comparison really fast. 

[`miniprot`](https://github.com/lh3/miniprot) spliced alignment of protein seq to genome seq (similar to genewise and exonerate). Not optimized for distant comparisons but it's doable.

## EMBL-EBI ProtVar database

https://www.ebi.ac.uk/ProtVar

They've done functional predictions on every possible human missense variant.

Many different types of variant predictor.
- view on structure
- score likelihood of having detrimental effect, e.g. AlphaFold missense score, CADD score, stability change


## Proteomics

[π-HuB: the proteomic navigator of the human body](https://www.nature.com/articles/s41586-024-08280-5)


# Electronic lab notebooks

LabArchive

# Not biology or coding

## Drawing programs

Some free alternatives to Adobe Illustrator:
- inkscape (Erick Matsen recommended)
- gimp
- paint.net
- AffinityDesigner Shirleen recommended, says it costs a $50 (one-time, not subscription) 



# Coding

## R
"escape excel" tool to get tab-delimited files ready to import into excel, without gene names getting mangled
https://github.com/pstew/escape_excel

## Python
Evolution python tools: http://etetoolkit.org/

## Markdown
To print a markdown file, use VScode. Right-click when the file is open to get a menu.
