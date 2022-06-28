#!/bin/bash
source /app/lmod/lmod/init/profile
module load SAMtools/1.11-GCC-10.2.0
module load BWA/0.7.17-GCC-10.2.0

(bwa mem -t 4 REFGENOME.fa mySample_R1.fq.gz mySample_R2.fq.gz | samtools view -@ 4 -Sb -  | samtools sort -@ 4 -O bam -T ${DELETE30}/malik_h/mySample_tmpDir > ./mySample.bwa.bam ) 2>> ./mySample.bwa.log.txt

samtools index ./mySample.bwa.bam

samtools flagstat ./mySample.bwa.bam > ./mySample.bwa.bam.flagstats

module purge
