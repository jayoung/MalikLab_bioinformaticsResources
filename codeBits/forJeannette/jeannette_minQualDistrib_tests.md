


xxxx test for jeannette project.  i want to see distribution of min qual each read

picard QualityScoreDistribution - not useful. shows all quals
```
cd ~/FH_fast_storage/bat_reproduction/data/fastq/artibeus_jamaicensis/dna/illumina/Orr_batch1/mapping/BWA_refseqAssembly
module purge
module load picard/2.25.0-Java-11
module load fhR/4.2.0-foss-2021b
java -Xmx10g -Djava.io. -Djava.io.tmpdir=${DELETE30}/malik_h -jar $EBROOTPICARD/picard.jar QualityScoreDistribution \
      -I Hch08.combined.testRegion100kb.bam \
      -O Hch08.combined.testRegion100kb.qual_score_dist.txt \
      -CHART Hch08.combined.testRegion100kb.qual_score_dist.pdf
module purge
```

fastp? https://github.com/OpenGene/fastp
no
```
module purge
module load fastp/0.20.0-GCC-8.3.0

fastp -i DMSLibraryPrep1_BBmapMerge.fixNames.fq.gz  -o DMSLibraryPrep1_BBmapMerge.fixNames.fastp.fq.gz -h DMSLibraryPrep1_BBmapMerge.fixNames.fastp.html -j DMSLibraryPrep1_BBmapMerge.fixNames.fastp.json

module purge
```

