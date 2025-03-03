# Genomics miscellaneous notes

## Definitions

"biclustering" - cluster genes and samples simultaneously, rather than independently (I think?)

## File formats

### sam/bam format

[Sam format specifications](http://samtools.github.io/hts-specs/SAMv1.pdf) and [optional fields](https://samtools.github.io/hts-specs/SAMtags.pdf)

## Algorithms

### Aligners

People are working on ways to align query sequences to much larger datasets of reference genomes, e.g. for metagenomics studies.

E.g. LexicMap (Zamin Iqbal lab) aligns queries >500bp long to millions of prokaryotic genomes much more quickly than BLAST

### Variant calling

#### Large variants (segmental variation, SVs)


- Delly and Sniffles (both use mapping to reference genome)
- SVarp for graph-based SV calling
- giggles for long-read based pangenome-graph-based re-genotyping

### Other

Clustering very large databases of sequences - [MMseqs2](https://github.com/soedinglab/MMseqs2) for protein or DNA. See [paper](https://www.nature.com/articles/nbt.3988). It's much faster than BLAST. Starts with kmer analysis, progresses to alignments.


## Projects/resources

xxxx I'm sure I have other lists like this elsewhere that I need to combine

There's a rapid release version of the Ensembl site that might have more than the main site.

### Human variation

[All Of Us](https://allofus.nih.gov) 
- NIH project
- 1m people living in US
- health focussed
- data: genomics, health metrics, surveys

### Expression

#### Databases/browsers for expression data

GTEX

TCGA

Human Protein Atlas



##### Kaessman lab

Kaessmann Lab's [EvoDevo mammalian organs](https://apps.kaessmannlab.org/alternative-splicing/) app. Uses their data from 8 tissues across development, 7 species (human, macaque, mouse, rat, rabbit, opossum, chicken).

Kaessman lab has other shiny apps:
- SpermEvol: "allows the interactive exploration of individual gene expression profiles of adult testis across ten mammalian species and a bird at the single-cell resolution.". Data from this paper [Murat, Mbengue et al. The molecular evolution of spermatogenesis across mammals (2021)](https://www.nature.com/articles/s41586-022-05547-7)
- SexBiasApp: "Evo-devo sex-biased mammalian genes" allows the interactive exploration of sex-biased gene expression profiles across organs, developmental stages and species.
- Ex2plorer: from the publication: [Wang, Z.Y., Leushkin, E. et al. Transcriptome and translatome co-evolution in mammals. Nature (2020)](https://www.nature.com/articles/s41586-020-2899-z)


#### Recount project 

The `recount` project aims to uniformly process all public RNA-seq data (for human only?)

`recount2` pipeline:
- published 2017 in [Nature Biotechnology](https://www.nature.com/articles/nbt.3838)
- mapping with [Rail-RNA](https://docs.rail.bio/tutorial/) ([published](https://academic.oup.com/bioinformatics/article/33/24/4033/2525684) in 2016) against hg38 generates
    - "cross-sample tables" (tsv)
    - bigwig coverage data
    - junction files (jx format)
- that gets wrangled to a `RangedSummarizedExperiment` R object
- more info in the [recount-contributions](https://github.com/leekgroup/recount-contributions) git repo


[`recount3`](https://rna.recount.bio):
- human and mouse
- ~ 10,000 studies for each species
- "raw sequencing data were processed with the Monorail system as described in the recount3 paper which created the coverage bigWig files and the recount-unified text files."
- published in [Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02533-6)
- Sophie Kogut told me that recount3 R package might not have all the functions that are in previous recount package versions.

There's an R package called recount - see [vignette](https://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html). There's a separate `recount3` bioconductor package

There's a [recount shiny app](https://jhubiostatistics.shinyapps.io/recount/). Not sure what it does, but I think it lets you find useful datasets and download gene expression.

Sophie also told me that Daniel Blanco-Melo has had advice from Leo (Leonardo Collado-Torres, recount author), and Monica from their lab has been using recount a bit. Also Jeff Leek was involved in this project (unclear whether he still is).

#### Single cell expression

[Human Common Cell Atlas](https://www.humancellatlas.org) is a repository of uniformly processed public single cell data