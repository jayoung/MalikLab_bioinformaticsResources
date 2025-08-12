require(ggplot2)
require(tidyr)
require(dplyr)
require(ggpubr)

## on Mac
masterH2Bdir <- "/Volumes/malik_h/user/jayoung/forOtherPeople/forPravrutha/H2B"
## on server
# masterH2Bdir <- "/fh/fast/malik_h/user/jayoung/forOtherPeople/forPravrutha/H2B" 

source(paste(masterH2Bdir, "H2B_functions.R", sep="/")) ## moved the color scheme into this file

#### old-school plots. one plot per gene, all samples in the dataset
plotRPKMsEachGene <- function(rpkmTable, plotFile, 
                              pageWidth=11, pageHeight=15,
                              textToStripFromSampleNames=NULL,
                              numPlotsHorizontal=2, numPlotsVertical=5,
                              minYaxisUpperLimit=5, # if no RPKM is >5, ylim[2] will be 5, 
                              # otherwise it talks max(rpkm)
                              pageTitles=NULL,
                              myOMA=c(5,0,5,0)
) {
  if (dim(rpkmTable)[1]>100) {
    stop("\n\n    You asked to make plots, but there are >100 genes - this is unwise, and I'm skipping it. Edit the script to override that choice\n\n")
  } 
  pdf(width=pageWidth,height=pageHeight, file=plotFile)
  par(mfrow=c(numPlotsVertical,numPlotsHorizontal), oma=myOMA)
  myNumPlotsPerPage <- numPlotsVertical * numPlotsHorizontal
  plotCount <- 0
  for (i in 1:dim(rpkmTable)[1]) {
    thisGene <- as.character(rpkmTable[i,"name"])
    cat ("    plotting RPKMs for",thisGene,"\n")
    thisGeneDat <- rpkmTable[i,grep("RPKM", colnames(rpkmTable))]
    
    sampleNames <- gsub("_RPKM","",colnames(thisGeneDat))
    if(!is.null(textToStripFromSampleNames)) {
      sampleNames <- gsub(textToStripFromSampleNames,"",sampleNames)
    }
    thisGeneDat <- as.numeric(thisGeneDat)
    names(thisGeneDat) <- sampleNames
    
    thisPlotRange <- range(thisGeneDat)
    thisPlotRange[1] <- 0
    if (max(thisGeneDat)<minYaxisUpperLimit) {thisPlotRange[2] <- minYaxisUpperLimit}
    barplot(thisGeneDat , ylim=thisPlotRange, main=thisGene, las=2, ylab="reads per million reads mapped")
    plotCount <- plotCount + 1
    if (plotCount == myNumPlotsPerPage) {
      title(main=pageTitles, outer=TRUE)
      plotCount <- 0
    }
  }
  if(!is.null(pageTitles)) { title(main=pageTitles, outer=TRUE) }
  dev.off()
} 


##### collectRPKMsForGGplot is a function to take selected samples, perhaps from >1 dataset, and get just those RPKMs in a nice format to plot
# sampleList is a list that looks like samplesOfInterest. First key is species, second key is dataset name, then we have a character vector of the samples to look at.
collectRPKMsForGGplot <- function(sampleList, species, metadata=rpkms_metadata, rpkmTables=rpkms) {
  if (!species %in% names(sampleList)) {
    stop("\n\nERROR - species ",species," is not in the sample list provided\n\n")
  }
  samples <- sampleList[[species]]
  
  ## get the relevant rpkm table names and the tables
  if (sum(!names(samples) %in% metadata[,"shortName"])>0) {
    missingSamples <- names(samples)[ which (!names(samples) %in% metadata[,"shortName"]) ]
    stop("\n\nERROR - there are RPKM tables that are not listed in the metadata table:",missingSamples,"\n\n")
  }
  datasetNames <- metadata[match(names(samples), metadata[,"shortName"]),"dataset"]
  names(datasetNames) <- names(samples)
  
  if (sum(!datasetNames %in% names(rpkmTables))>0) {
    missingDatasets <- datasetNames[ which (!datasetNames %in% names(rpkmTables)) ]
    stop("\n\nERROR - there are datasets that are not present in the rpkm tables I read in ",missingDatasets,"\n\n")
  }

  ## get RPKM data for the selected samples 
  rpkmTablesSelectedSamples <- lapply(names(samples), function(x) {
    myColNames <- paste(samples[[x]], "_RPKM", sep="")
    myTable <- rpkmTables[[ datasetNames[x] ]]
    if(sum(!myColNames %in% colnames(myTable))>0) {
      missingColnames <- myColNames[ which(!myColNames %in% colnames(myTable)) ]
      missingColnames <- paste(missingColnames, collapse=" ")
      stop("\n\nERROR - some colnames specified for dataset ",x," are not present in the corresponding RPKM table - these are the missing column names:\n",missingColnames,"\n\n")
    }
    selectedDat <- myTable[,myColNames, drop=FALSE]  ## the drop=FALSE ensure I get a data.frame back (with a colname) even if there is only a single sample
    
    ## xx if there's only one sample it returns it as a vector rather than a data.frame
    return(selectedDat)
  })
  #names(rpkmTablesSelectedSamples) <- names(samples)
  ## combine RPKMs in a single table
  if (length(samples)>1) {
    newDat <- do.call("cbind", rpkmTablesSelectedSamples)
  } else {
    newDat <- rpkmTablesSelectedSamples[[1]]
  }
  rownames(newDat) <- rpkmTables[[ datasetNames[1] ]] [,"name"]
  return(newDat)
}

#### convertToLongFormat - use this on a table where rows are genes, rownames are gene names, columns are ALL RPKMs. Converts it to a long format tibble
convertToLongFormat <- function(rpkmOnlyTable) {

  rpkmOnlyTable_long <- pivot_longer(as.data.frame(t(rpkmOnlyTable)), 
                                     cols=everything(), 
                                     names_to = "Gene", values_to = "RPKM")
  rpkmOnlyTable_long[,"sample"] <- rep(colnames(rpkmOnlyTable), each=dim(rpkmOnlyTable)[1])
  rpkmOnlyTable_long <- as.data.frame(rpkmOnlyTable_long)
  return(rpkmOnlyTable_long)                      
}
