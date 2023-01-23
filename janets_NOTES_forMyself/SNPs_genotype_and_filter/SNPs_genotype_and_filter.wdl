version 1.0

## SNPregionBEDfiles comes from 
# GCF_014825515.1_WHU_Ajam_v2_genomic.bedList.small.json (2 chunks)
# GCF_014825515.1_WHU_Ajam_v2_genomic.bedList.full.json (121 chunks)

## VCFfiles comes from 
# VCFfileNames.twoSamplesTwoChunks.json 

## commonInputs comes from SNPs_genotype_and_filter.commonInputs.json

## define a structure for variables in the common inputs json
struct commonInputsStruct {
  File refGenomePicard_fasta
  File refGenomePicard_fai
  File refGenomePicard_dict
  File script_VCFgetSomeCtgs
  File refGenome_bigCtgNamesFile
}


workflow SNPs_genotype_and_filter {
  input {
    Array[String] allSamples
    Array[String] allChunks
    Map[String,Map[String,Array[File]]] VCFfiles 
    Array[Pair[String,File]] SNPregionBEDfiles
    #### commonInputs (reference genome, etc)
    commonInputsStruct commonInputs
  }

  ## scatter over each chunk to run GenomicsDBImport and GenotypeGVCFs (for each chunk, we use the chunk files from all samples)
  scatter (thisRegion in SNPregionBEDfiles) {
    String regionName = thisRegion.left
    File regionFile = thisRegion.right

    call GenomicsDBImport {
      input:
        chunkName = regionName,
        bedFile = regionFile,
        vcfs = VCFfiles[regionName]['vcfFiles'],
        vcfIndices = VCFfiles[regionName]['vcfIndexFiles']
    } 
    call GenotypeGVCFs {
      input:
        gendb = GenomicsDBImport.genotypeDB,
        chunkName = regionName,
        bedFile = regionFile,
        refGenome_fa = commonInputs.refGenomePicard_fasta,
        refGenome_fai = commonInputs.refGenomePicard_fai,
        refGenome_dict = commonInputs.refGenomePicard_dict
    }
  } # end of eachChunk scatter

  ## merge vcfs from all genomic chunks
  call mergeVCFs {
    input:
      genotypedVCFs = GenotypeGVCFs.genotypeVCF
  }

  ## split vcf into two files: SNPs and indels
  call split_SNPs_indels {
    input:
      VCF = mergeVCFs.mergedVcf
  }

  ## filter indels
  call filterIndels {
    input: 
      VCF = split_SNPs_indels.indelOut
  }

  ## filter SNPs
  call filterSNPs {
    input: 
      VCF = split_SNPs_indels.snpOut
  }
  ## get just biallelic SNPs
  call getBiallelicSNPs {
    input: 
      VCF = filterSNPs.filteredVCF
  }

  ## get only SNPs on large contigs
  call vcfGetOnlyBigContigs {
    input:
      getCtgsScript = commonInputs.script_VCFgetSomeCtgs,
      ctgNames = commonInputs.refGenome_bigCtgNamesFile,
      VCF = getBiallelicSNPs.biallelicVCF_shortHeader
  }
  # split vcf into one per ctg and run bcftools stats


  # R processing on individual contig vcf files

  # get sites that vary only within the Orr colony


  ### we're done!
  # xx make sure this is output of whatever is the last task
  call send_finished_message {
    input:
      fileNames = vcfGetOnlyBigContigs.bigCtgsVCFidx
  }

  ## final outputs to emit from the workflow
  # xx need to add more outputs here
  output {
    #Array[File] genotypeDB = GenomicsDBImport.genotypeDB
    #Array[File] genotypedVCF = GenotypeGVCFs.genotypeVCF
    File mergedVcf = mergeVCFs.mergedVcf
    File mergedVcfIdx = mergeVCFs.mergedVcfIdx
    File mergedVcfLog = mergeVCFs.mergedVcfLog
  }

}



######### task definitions

task GenomicsDBImport {
  input {
    String chunkName
    File bedFile
    Array[File] vcfs
    Array[File] vcfIndices
  }
  command <<<
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      --genomicsdb-workspace-path ~{chunkName}.genotyping.db \
      --tmp-dir /fh/scratch/delete10/malik_h/jayoung \
      --merge-contigs-into-num-partitions 2 \
      -L ~{bedFile} \
      -V ~{sep=' -V ' vcfs} >> ~{chunkName}.GenomicsDBImport.log.txt 2>&1
  >>>
  runtime {
    modules: "GATK/4.2.6.1-GCCcore-11.2.0"
    memory: "5GB"
  }
  output {
    File genotypeDB = "~{chunkName}.genotyping.db"
  }
}

task GenotypeGVCFs {
  input {
    File gendb
    String chunkName
    File bedFile
    File refGenome_fa
    File refGenome_fai
    File refGenome_dict
  }
  command <<<
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R ~{refGenome_fa} \
      -V gendb://~{gendb} \
      --tmp-dir /fh/scratch/delete10/malik_h/jayoung \
      -O ~{chunkName}.genotypes.vcf.gz >> ~{chunkName}.genotypes.log.txt 2>&1
  >>>
  runtime {
    modules: "GATK/4.2.6.1-GCCcore-11.2.0"
    memory: "5GB"
  }
  output {
    File genotypeVCF = "~{chunkName}.genotypes.vcf.gz"
    File genotypeLog = "~{chunkName}.genotypes.log.txt"
    #File genotypeVCFidx = "~{chunkName}.genotypes.vcf.gz.tbi"
  }
}

task mergeVCFs {
  input {
    Array[File] genotypedVCFs
  }
  command <<<
    java -Xmx10g -Djava.io. \
      -Djava.io.tmpdir=/fh/scratch/delete10/malik_h/jayoung \
      -jar $EBROOTPICARD/picard.jar MergeVcfs \
      -I ~{sep=' -I ' genotypedVCFs} \
      -O allChunks_combined.vcf.gz >> allChunks_combined.MergeVcfs.log.txt 2>&1
  >>>
  runtime {
    modules: "picard/2.25.0-Java-11"
    memory: "12GB"
  }
  output {
    File mergedVcf = "allChunks_combined.vcf.gz"
    File mergedVcfIdx = "allChunks_combined.vcf.gz.tbi"
    File mergedVcfLog = "allChunks_combined.MergeVcfs.log.txt"
  }
}

task split_SNPs_indels {
  input {
    File VCF
  }
  command <<<
    java -Xmx10g -Djava.io. \
      -Djava.io.tmpdir=/fh/scratch/delete10/malik_h/jayoung \
      -jar $EBROOTPICARD/picard.jar SplitVcfs \
      I=~{VCF} \
      SNP_OUTPUT=allChunks_combined.snp.vcf \
      INDEL_OUTPUT=allChunks_combined.indel.vcf \
      STRICT=false
  >>>
  runtime {
    modules: "picard/2.25.0-Java-11"
    memory: "12GB"
  }
  output {
    File snpOut = "allChunks_combined.snp.vcf"
    File indelOut = "allChunks_combined.indel.vcf"
  }
}


task filterSNPs {
  input {
    File VCF
  }
  command <<<
    gatk VariantFiltration \
        -V ~{VCF} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O allChunks_combined.snps.addFilt.vcf.gz
    # convert to table, for R (including variants that failed the filtering):
    gatk VariantsToTable \
        -V allChunks_combined.snps.addFilt.vcf.gz \
        -F CHROM -F POS -F TYPE -F QUAL -F FILTER -F QD -F FS \
        -F SOR -F MQ -F MQRankSum -F ReadPosRankSumTest \
        --show-filtered TRUE \
        -O allChunks_combined.snps.addFilt.vcf.table.tsv 
  >>>
  runtime {
    modules: "GATK/4.2.6.1-GCCcore-11.2.0"
  }
  output {
    File filteredVCF = "allChunks_combined.snps.addFilt.vcf.gz"
    File filteredVCFtsv = "allChunks_combined.snps.addFilt.vcf.table.tsv"
  }
}

task getBiallelicSNPs {
  input {
    File VCF
  }
  command <<<
    # get only biallelic SNPs
    bcftools view -f 'PASS,.' --max-alleles 2 --output-type z ~{VCF} > allChunks_combined.snps.addFilt.biallelicFilt.vcf.gz
    bcftools index allChunks_combined.snps.addFilt.biallelicFilt.vcf.gz

    # get version without header
    zcat allChunks_combined.snps.addFilt.biallelicFilt.vcf.gz | grep -v '^##contig' | bgzip > allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.vcf.gz
    bcftools index allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.vcf.gz
  >>>
  runtime {
    modules: "BCFtools/1.14-GCC-11.2.0"
  }
  output {
    File biallelicVCF_shortHeader = "allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.vcf.gz"
    File biallelicVCFidx_shortHeader = "allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.vcf.gz.csi"
  }
}

## filter the no-header-biallelic SNPs vcf for only variants on contigs of at least 100kb size
task vcfGetOnlyBigContigs {
  input {
    File getCtgsScript
    File ctgNames
    File VCF
  }
  command <<<
    ~{getCtgsScript} --outdir=. --contigs=~{ctgNames} ~{VCF}
  >>>
  runtime {
    modules: "BCFtools/1.14-GCC-11.2.0"
  }
  output {
    File bigCtgsVCF = "allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.someCtgs.vcf.gz"
    File bigCtgsVCFidx = "allChunks_combined.snps.addFilt.biallelicFilt.noCtgHdr.someCtgs.vcf.gz.csi"
  }
}

## split the vcf, making one for each contig
# task splitVCFeachContig {

# }


task filterIndels {
  input {
    File VCF
  }
  command <<<
    gatk VariantFiltration \
        -V ~{VCF} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O allChunks_combined.indels.addFilt.vcf.gz >> allChunks_combined.indels.addFilt.log.txt 2>&1 
    # convert to table, for R (including variants that failed the filtering):
    gatk VariantsToTable \
        -V allChunks_combined.indels.addFilt.vcf.gz \
        -F CHROM -F POS -F TYPE -F QUAL -F FILTER -F QD -F FS \
        -F SOR -F MQ -F MQRankSum -F ReadPosRankSumTest \
        --show-filtered TRUE \
        -O allChunks_combined.indels.addFilt.vcf.table.tsv 
  >>>
  runtime {
    modules: "GATK/4.2.6.1-GCCcore-11.2.0"
    memory: "12GB"
  }
  output {
    File filteredVCF = "allChunks_combined.indels.addFilt.vcf.gz"
    File filteredVCFtsv = "allChunks_combined.indels.addFilt.vcf.table.tsv"
  }
}


task send_finished_message {
  input {
    File fileNames
  }
  # File firstFile = fileNames[0]
  command <<<
    echo "SNPs_genotype_and_filter workflow has finished.\nFirst of the final output files is called:\n~{fileNames}." | mail -s "SNPs_genotype_and_filter workflow has finished" jayoung@fredhutch.org
  >>>
  output {
    File stdout = stdout()
  }
}

