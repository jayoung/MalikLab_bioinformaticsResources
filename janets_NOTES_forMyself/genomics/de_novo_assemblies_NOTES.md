# de novo genome (and transcriptome) assemblies

## Definitions

"contig" versus "unitig" - "Unitigs, formed by kmers, represent the sequences without branch in assembly graph. Several connected unitigs form a contig and one unitig can be present in multiple contigs."

"MAG" metagenome-assembled genome - microbial genomes reconstructed from metagenome data.

"PVG" pangenome variation graphs
- can map reads to PVGs (using special tools) to mitigate reference bias

"DBG" de Bruijn graph

## Assemblies

Algorithms:
- metaMDBG - [github](https://github.com/GaetanBenoitDev/metaMDBG) and [paper](https://www.nature.com/articles/s41587-023-01983-6) a metagenomics assembler for PacBio HiFi reads (and I think ONT reads)

Often for assembly there is a stage of graph trimming (loose ends) and error correction

## Pangenome and variation graphs

Full variation graph (every variant, including SNPs) might be bigger than a pangenome graph, which more likely only branches for larger indels and structural variants.

Someone I met at the 2024 Genomic Informatics meeting Pierre Peterlongo told me to check out recent pangenome PhD by Sandra Romain

Papers
- [PanGenome Graphs](https://pmc.ncbi.nlm.nih.gov/articles/PMC8006571/) review article, Eizenga et al 2021 Annu Rev Genomics Hum Genet


Resources with many links
- https://pangenome.github.io
- [ALPACA](https://alpaca-itn.eu) - ALgorithms for PAngenome Computational Analysis

Algorithms:
- [minigraph toolkit](https://pubmed.ncbi.nlm.nih.gov/33066802/)
- minigraph cactus - [github](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) and [paper](https://www.nature.com/articles/s41587-023-01793-w). creates pangenomes directly from whole-genome alignments. Used on human. For more diverged haplotypes they recommend [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md)
- PGGB pangenome graph
- Panaroo (gene-based? gene family based?)
- [Metagraph](https://github.com/ratschlab/metagraph) - a tool for scalable construction of annotated genome graphs and sequence-to-graph alignment.
- [ODGI tools](https://academic.oup.com/bioinformatics/article/38/13/3319/6585331). Optimized Dynamic Genome/Graph Implementation. Suite of tools for graph analysis.

- graphtools (sounds a bit like samtools for graphs)

- SVarp for graph-based SV calling


Visualizing assembly graphs
- Bandage (and perhaps can give summary stats?)
- SequenceTubeMap (for pangenome variation graphs)
- WARAGRAPH


See also algorithms shown [here](https://github.com/jayoung/thoughts/blob/main/notes/conferences/2024_GenomeInformatics/2024_GenomeInformatics_NOTES.md#tim-downing-poster---viral-pangenomes)


## using pangenome graphs

PanGenie - match short reads to a pangenome graph for genotyping known variants