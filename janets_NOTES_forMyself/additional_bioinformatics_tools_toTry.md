# Getting species trees

[TimeTree](http://timetree.org)

[phyloT](https://phylot.biobyte.de) for getting accepted species trees


# Phylogenies


## IQ-Tree

I see this used in several large-scale studies

## Treebest

[Treebest](http://treesoft.sourceforge.net/treebest.shtml) - used in Ensembl pipeline. Combines multiple tree types (nucleotide, protein, dN, dS, species tree), tries to reduce number of mismatches to species tree, infers duplication, deletion, etc.

# Phylogenetic tree display/plotting

interactive Tree Of Life [iTOL](https://itol.embl.de) has some nice-looking plotting tools

[ggtree](https://bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/ggtree.html) for R

# Transposons

TEtools - [paper](https://www.ncbi.nlm.nih.gov/pubmed/27924026) and [repo](https://github.com/l-modolo/TEtools)

TEanalysis

RepEnrich2


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

# Datasets/resources

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
