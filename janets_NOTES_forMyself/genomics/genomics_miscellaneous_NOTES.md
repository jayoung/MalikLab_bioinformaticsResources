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

#### Single cell expression

[Human Common Cell Atlas](https://www.humancellatlas.org) is a repository of uniformly processed public single cell data