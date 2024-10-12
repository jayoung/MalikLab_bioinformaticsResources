# search tools

NCBI's Edirect toolkit to help search and extract data from their databases:
https://www.ncbi.nlm.nih.gov/books/NBK179288/
https://ncbi-hackathons.github.io/EDirectCookbook/

NCBI's SRAtoolkit

# Getting species trees

[TimeTree](http://timetree.org)

[phyloT](https://phylot.biobyte.de) for getting accepted species trees

NCBI taxonomy

Open Tree of Life
https://tree.opentreeoflife.org/opentree/argus/opentree9.1@ott93302

Lifemap displays of NCBI taxonomy, or the open tree of life:
http://lifemap.univ-lyon1.fr/

# Phylogenies

IQ-Tree - I see this used in several large-scale studies

[Treebest](http://treesoft.sourceforge.net/treebest.shtml) - used in Ensembl pipeline. Combines multiple tree types (nucleotide, protein, dN, dS, species tree), tries to reduce number of mismatches to species tree, infers duplication, deletion, etc.

NOTUNG for tree reconciliation

CAFE for tree reconciliation (David Liberles says may be too simplistic)

David Liberles: expectation maximization for missing data

SATe - simultaneous alignment and tree estimation - Liu, Raghavan, Nelesen, Linder and Warnow. Science, 2009

Fast Statistical Alignment
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000392
http://orangutan.math.berkeley.edu/fsa/, http://fsa.sourceforge.net/.

PLATO for detecting variation in evolutionary rates along alignment

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


# user-friendly online toolsets

Galaxy
- Main public instance
- Hutch instance
- various other instances

EDGE - seems to be mostly for bacterial / metagenomics analyses?

# protein sequence analysis

for structure, Antoine used UCSF Chimera (Pettersen, et al. 2004).

# Deep sequencing

GC normalization: John Huddleston advised using mrFAST or mrsFAST for read mapping and then mrcanavar for GC normalization and copy number estimation. But that was a while ago (<2017?)

## RNA-seq
Ching-Ho likes HISAT2 (new version of tophat) for RNA-seq., partly because it has fewer options than STAR, but he also thinks STAR might do better with multi-mapping genes.  HISAT2 may be able to directly take SRA accessions and map without downloading reads.

Ching-Ho also notes that some RNA-seq datasets have high rRNA contamination. This could be affecting my RPKMs quite a lot.


# Coding

## R
"escape excel" tool to get tab-delimited files ready to import into excel, without gene names getting mangled
https://github.com/pstew/escape_excel

## Python
Evolution python tools: http://etetoolkit.org/

## Markdown
To print a markdown file, use VScode. Right-click when the file is open to get a menu.

# Datasets/resources

## Genome browsers

UCSC
for custom tracks, hubs, etc, we can use a directory on fast that is accessible from beyond the Hutch firewall: /fh/fast/malik_h/pub/
see the notes file in that directory for more detail on how to host a hub (malik_pubDir_NOTES_JY.txt)

Ensembl

## species-specific databases

Flybase, Wormbase, SGD

## Orthology / conservation databases and data sources

OrthoDB (there is an R package for OrthoDB). Provides estimate of evolutionary rates

[OrthoMaM](https://academic.oup.com/nar/article/52/D1/D529/7318103) - 15,868 alignments from 190 mammalian genomes, should be 1-to-1 othologs. Nice alignment filtering: should be high quality? Should be in-frame.

Zoonomia project: >240 mammals analyzed. Precomputed genome-wide mammalian conservation scores, and codon alignments are here https://zoonomiaproject.org/the-data/ 

12-species alignments (10 mammals+chicken+tortoise) (nucleotide and amino acid) for 7445 genes from Shuler & Hagai paper 2022 in teh associated [github repo](https://github.com/galshuler/Evolution_of_host-virus_protein_protein_interactions). I downloaded and unpacked in `~/FH_fast_storage/paml_screen/Adrian_Pelin_vaccinia/references/2022-Shuler/git_repo/Evolution_of_host-virus_protein_protein_interactions`. See notes in `~/FH_fast_storage/paml_screen/Adrian_Pelin_vaccinia/Adrian_Pelin_vaccinia_NOTES.md`

TreeFam

COG: NCBI's database of clusters of orthologous genes.  [COG website](https://www.ncbi.nlm.nih.gov/research/cog) allows you to send selected protein sequences from COG to COBALT alignment generator, to tweak alignment parameters if needed, and then to download an alignment.

PhyloPat
Seems to be out-of-date and unmaintained
June 2014: http://www.cmbi.ru.nl/cdd/phylopat/52/
Based on Ensembl 52 (quite old:  Ensembl now up to v75)

[Consurf-db](https://consurfdb.tau.ac.il/overview.php) 
- per-residue rates for proteins that have structure in PDB. 
- Uses up to 300 homologues, found in UniRef90, with some seq identity/coverage filters 
- evolutionary rates for each site calculated based on tree+alignment using Rate4Site algorithm

[orthogene](https://bioconductor.org/packages/release/bioc/vignettes/orthogene/inst/doc/orthogene.html) R package


phyloP tracks from UCSC genome browser



## Methods for conservation / evolutionary rates

### nucleotide

dN/dS type methods

genome-based methods (phyloP etc)

[SLR](https://www.ebi.ac.uk/research/goldman/software/slr/)

Selecton, 2007 - [publication](https://pubmed.ncbi.nlm.nih.gov/17586822/) and [website](http://selecton.bioinfo.tau.ac.il). Website not responding on 10/11/2024. I don't see a way to download the algorithm. Looking at the paper, it looks like it uses various evolutionary models (e.g. codeml M7/M8 etc, or MEC, which takes protein models like JTT and converts them to codon models). Does some kind of hypothesis testing, chooses the best-fitting model.  Also uses a Bayesian method to estimate sitewise dN/dS - uses 8 discrete rate categories and does something with posteriors to get the estimate.

### protein

[Evorator](https://evorator.tau.ac.il) - website not responding on 10/4/24

[Rate4Site](https://www.tau.ac.il/~itaymay/cp/rate4site.html) (ref [Pupko et al, 2002](https://pubmed.ncbi.nlm.nih.gov/12169533/)) - given a multiple sequence alignment, it can generate a tree (or you can provide one) and it calculates the relative evolutionary rate at each site.

## Transcriptomics

Multispecies developmental timecourse [browser](https://apps.kaessmannlab.org/evodevoapp/).
- human, rhesus macaque, mouse, rat, rabbit, opossum and chicken
- seven organs (cerebrum, cerebellum, heart, kidney, liver, ovary and testis)
- developmental time points from early organogenesis to adulthood

## Single cell data

Brotman Baty's [Descartes](https://descartes.brotmanbaty.org) - gateway to single cell atlas for several species

## Variation 

Human: 1000 genomes [selection browser](http://hsb.upf.edu/)

Drosophila: [popfly](http://popfly.uab.cat/) (has data to allow M-K tests)



# Electronic lab notebooks

LabArchive

# Not biology

## Drawing programs

Some free alternatives to Adobe Illustrator:
- inkscape (Erick Matsen recommended)
- gimp
- paint.net
- AffinityDesigner Shirleen recommended, says it costs a $50 (one-time, not subscription) 
