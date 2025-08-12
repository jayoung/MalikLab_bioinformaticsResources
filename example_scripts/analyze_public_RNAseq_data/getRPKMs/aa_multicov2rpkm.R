##### takes 
# (1) first arg = multicov file 
# (2) other args = flagstats files for each sample in the multicov file 
#    calculates RPKM for each gene, using total number of mapped reads as denominator
#    calculates percentiles for those RPKMs, only considering RPKMs that are above the threshold
#    for hg19 and mm10, I translate RefSeq IDs to gene names using tables I downloaded from UCSC

## xxx I could also choose just a single transcript for each gene before taking RPKMs.
## xxx I could also use alias files to add SRA IDs for each sample names

### usage: 
# Rscript /home/jayoung/Rscripts/multicov2rpkm.R muticovFilename.txt flagstatsFiles.txt >& multicovFilename.Rout

RPKMthreshold <- 1

getPercentiles <- FALSE
## In some cases (e.g. files with not many genes) I don't want percentiles, and they might create errors (e.g. if I have one gene, and some tissues in which it has sub-threshold RPKM - then there are no remaining values from which to calculate percentiles)

makePlots <- FALSE
# if yes, will make a plot for all genes, using all samples

##########################

source ("/fh/fast/malik_h/user/jayoung/Rscripts/JY_bigseq_functions.R")

if(getPercentiles) {
    cat("\nFor percentiles, I ignored genes with RPKM <",RPKMthreshold,"\n\n")
} else {
    cat("\nNot getting percentiles - flag is set within script\n\n")
}

### get file names from command line
allargs <- commandArgs(TRUE)
countsFile <- allargs[1]
flagstatsFiles <- allargs[2:length(allargs)]

# temp 
# countsFile <- "Plat_genes.blat_ornAna2.besthits.psl.bed.counts.platypus_Brawand.split.multicov"
# flagstatsFiles <-list.files("/fh/fast/malik_h/grp/public_databases/NCBI/SRA/data/mammalian_expression_profiles/platypus/platypus_Brawand/STAR_ornAna2_maxMultiHits1_renamed", pattern="flagstats", full.names=TRUE)


### read flagstats files
cat("reading flagstats files\n")
allFlagstats <- combineFlagstatsFiles( flagstatsFiles )
colnames(allFlagstats) <- gsub(".flagstats$","", colnames(allFlagstats))
colnames(allFlagstats) <- gsub("_tophat.bam$","", colnames(allFlagstats))
colnames(allFlagstats) <- gsub("_STAR.bam$","", colnames(allFlagstats))
colnames(allFlagstats) <- gsub(".merged.bam$","", colnames(allFlagstats))
if (grepl("/", colnames(allFlagstats)[1])) {
    colnames(allFlagstats) <- sapply( strsplit(colnames(allFlagstats), "/"), tail, 1)
}

### read counts files
cat("reading counts file",countsFile,"\n")
# the check.names=FALSE is to preserve hyphens in column names
counts <- read.delim(countsFile, check.names=FALSE)
colnames(counts) <- gsub("_tophat.bam$","", colnames(counts))
colnames(counts) <- gsub("_STAR.bam$","", colnames(counts))
sampleNames <- colnames(counts)[13:dim(counts)[2]]

### check that we have flagstats for all samples in the counts file
if (sum(!sampleNames %in% colnames(allFlagstats)) > 0) {
    missingSamples <- setdiff(sampleNames, colnames(allFlagstats))
    missingSamples <- paste(missingSamples, collapse="\n")
    cat("TERMINATING - there are samples that aren't represented by a flagstats file:\n",missingSamples,"\n")
    quit(save="no")
}

### get alternative gene names (maybe)

### set up a function to add gene name (rather than transcript name) to counts table so I can later take just one isoform per gene
# function that uses an ensembl transcript-to-geneName translation table to add gene names on to a data.frame. transcript IDs should be in "name" column.
addGeneID <- function(myDF, myGeneIDs, stripOffAfterDot=TRUE) {
    ## strip .1 .2 off name before looking it up:
    if (stripOffAfterDot) {
        tryTheseNames <- sapply( strsplit( as.character(myDF[,"name"]), "\\." ), "[[", 1)
    } else {
        tryTheseNames <- as.character(myDF[,"name"])
    }
    if (sum( tryTheseNames %in% rownames(myGeneIDs))==0) {
        stop("    ERROR - gene names in the counts table don't seem to be in the IDs table\n\n")
    }

    temp <- sum( tryTheseNames %in% rownames(myGeneIDs))
    cat("    found",temp,"of",length(tryTheseNames),"IDs in the lookup table\n")
    
    newNames <- as.character(myGeneIDs[ match(tryTheseNames, rownames(myGeneIDs)), "geneSymbol"])
    myDF[,"geneName"] <- newNames
    myDF[which(is.na(newNames)),"geneName"] <- as.character(myDF[which(is.na(newNames)),"name"])
    myDF
}
geneNameTranslationFile <- NULL
geneNameTranslationType <- NULL
if ( grepl("hg19",countsFile) ) {
    geneNameTranslationFile <- "/fh/fast/malik_h/grp/public_databases/UCSC/Feb2009/other_tracks/kgXref_hg19_2017_Nov10.txt"
    geneNameTranslationType <- "refseq"
}
if ( grepl("mm10",countsFile) ) {
    geneNameTranslationFile <- "/fh/fast/malik_h/grp/public_databases/UCSC/mouse_Dec2011/other_tracks/kgXref_mm10_2017_Nov29.txt"
    geneNameTranslationType <- "refseq"
}
if (!is.null(geneNameTranslationFile)) {
    if (geneNameTranslationType == "refseq") {
        if (!file.exists(geneNameTranslationFile)) {
            cat("\n\n    TERMINATING - file",geneNameTranslationFile,"does not exist\n\n")
            quit(save="no")
        }
        geneIDs <- read.delim(geneNameTranslationFile)
        geneIDs <- geneIDs[ which(geneIDs[,"refseq"]!=""),]
        geneIDs <- unique(geneIDs[,c("refseq","geneSymbol")])
        geneIDs[,"refseq"] <- as.character(geneIDs[,"refseq"])
        geneIDs[,"geneSymbol"] <- as.character(geneIDs[,"geneSymbol"])
        ### some refseqs really are associated with >1 gene symbol. I will arbitrarily pick one - I don't care which 
        geneIDs <- geneIDs[which(!duplicated(geneIDs[,"refseq"])),] 
        rownames(geneIDs) <- geneIDs[,"refseq"]
        counts <- addGeneID(counts, geneIDs, stripOffAfterDot=FALSE)
    } else {
        cat("\n\n    TERMINATING - not sure how to deal with name translation file of this type\n\n")
        quit(save="no")
    }
}



### start getting RPKMs
# keep the first 12 or 13 (constant) columns from the bed12/counts file
rpkms <- counts[,1:12]
if ("geneName" %in% colnames(counts)) {
    rpkms[,"geneName"] <- counts[,"geneName"]
}

# get total gene size, add that as a column
getTotalBlockSizes <- function(myDF) {
    blockSizes <- as.character(myDF[,"blockSizes"])
    myDF[,"totalBlockSize"] <- sapply( strsplit(blockSizes, ","), function(x) {sum(as.integer(x))})
    myDF
}
rpkms <- getTotalBlockSizes(rpkms)

# add the gene counts
rpkms[,sampleNames] <- counts[,sampleNames]

# for each sample, get RPKM and percentiles
tempPercentileTable <- data.frame(row.names=1:dim(counts)[1])
geneSizesKB <- rpkms[,"totalBlockSize"]/1000
for (tempSample in sampleNames) {
    tempNumMappedReads <- allFlagstats["mapped",tempSample]/1000000
    tempRPKMs <- (counts[,tempSample] / geneSizesKB) / tempNumMappedReads
    
    ## add read counts to table
    tempCountColname <- paste(tempSample, "_reads", sep="")
    colnames(rpkms)[which(colnames(rpkms)==tempSample)] <- tempCountColname
    ## add RPKM value to table
    tempRPKMcolname <- paste(tempSample, "_RPKM", sep="")
    rpkms[,tempRPKMcolname] <- tempRPKMs
    
    rm(tempNumMappedReads, tempRPKMs)
}

rpkms <- cbind(rpkms, tempPercentileTable)

rm(tempSample, tempPercentileTable)

#### write output
outFile <- paste(countsFile, ".RPKMs.csv", sep="")
write.csv(rpkms, outFile, row.names=FALSE)

if (makePlots) {
    if (dim(rpkms)[1]>100) {
        cat("    You asked to make plots, but there are >100 genes - this is unwise, and I'm skipping it. Edit the script to override that choice\n")
    } else {
        dir <- getwd()
        pageTitles <- paste(countsFile,"\n",dir)
        plotFile <- paste(countsFile,".rpkmPlots.pdf", sep="")
        pdf(width=11,height=15, file=plotFile)
        par(mfrow=c(5,2), oma=c(5,0,2,0))
        myNumPlotsPerPage <- 10
        plotCount <- 0
        for (i in 1:dim(rpkms)[1]) {
            thisGene <- as.character(rpkms[i,"name"])
            cat ("plotting RPKMs for",thisGene,"\n")
            thisGeneDat <- rpkms[i,grep("RPKM", colnames(rpkms))]
            
            sampleNames <- gsub("_RPKM","",colnames(thisGeneDat))
            thisGeneDat <- as.numeric(thisGeneDat)
            names(thisGeneDat) <- sampleNames
            
            thisPlotRange <- range(thisGeneDat)
            thisPlotRange[1] <- 0
            #if (max(thisGeneDat)==0) {thisPlotRange[2] <- 1}
            if (max(thisGeneDat)<5) {thisPlotRange[2] <- 5}
            barplot(thisGeneDat , ylim=thisPlotRange, main=thisGene, las=2, ylab="reads per million reads mapped")
            plotCount <- plotCount + 1
            if (plotCount == myNumPlotsPerPage) {
                title(main=pageTitles, outer=TRUE)
                plotCount <- 0
            }
        }
        title(main=pageTitles, outer=TRUE)
        dev.off()
    } 
}
    
