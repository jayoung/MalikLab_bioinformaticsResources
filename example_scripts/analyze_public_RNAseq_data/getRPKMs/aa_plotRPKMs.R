rm(list=ls())

## on server:
wd <- "/Volumes/malik_h/user/jayoung/forOtherPeople/forPravrutha/H2B/geneSeqs/allGenes_2021_Mar/bySpecies/blatToGenomes/countReads/STAR/"
setwd(wd)
workingDir <- gsub("/Volumes/malik_h/user/jayoung/forOtherPeople/forPravrutha/","",wd)
RPKM_csv_dir <- "aa_RPKMs"
plotsDirRaw <- "aa_plots/raw"
plotsDirPretty <- "aa_plots/pretty"

### load functions
source("aa_plotRPKMs_functions.R")


##### some lists to help me plot only certain genes/samples:
# selected tissues (brain, liver, kidney, heart, testis, ovaries and embryos)  and 
# selected species (Chicken Opossum, Dog, Pig, Mouse, Rhesus, Human)
# selected genes
genesOfInterest <- list()
genesOfInterest[["dog"]] <- c("H2B.L4_Dog", "H2BMW_Dog-a", "H2Bopp8C_Dog", "SubH2B_Dog", "TSH2B_Dog")
genesOfInterest[["pig"]] <- c("H2B.L4_Pig", "H2Bopp8C_Pig","SubH2B_Pig","TSH2B_Pig-a" )
genesOfInterest[["mouse"]] <- c("H2BMW_Mouse", "SubH2B_Mouse", "TSH2B_Mouse")
genesOfInterest[["human"]] <- c("H2B.L4_Hum", "H2BMW_Hum-H2BM", "H2BMW_Hum-FWT", "H2Bopp8C_Hum", "SubH2B_Hum-pseudo", "TSH2B_Hum")
genesOfInterest[["chicken"]] <- c("H2B.L4_Chick")
genesOfInterest[["opossum"]] <- c("H2Bopp8C_Oppos-Hit2", "SubH2B_Oppos")
genesOfInterest[["macaque"]] <- c("H2Bopp8C_Macaque", "SubH2B_Macaque")

## one list per plot we will make. initialize now but fill with values later
samplesOfInterest <- list()
samplesOfInterest[["dog"]] <- list()
samplesOfInterest[["chicken"]] <- list()
samplesOfInterest[["opossum"]] <- list()
samplesOfInterest[["macaque"]] <- list()

samplesOfInterest[["pig"]] <- list()
samplesOfInterest[["pig_bodyMap_batch1"]] <- list()
samplesOfInterest[["pig_oocyteEarlyEmbryo"]] <- list()

samplesOfInterest[["mouse"]] <- list()
samplesOfInterest[["mouseXueOocyteEarlyEmbryo"]] <- list()

samplesOfInterest[["human"]] <- list()
samplesOfInterest[["humanSpermatogenesisStudy1"]] <- list()
samplesOfInterest[["humanSpermatogenesisStudy2"]] <- list()
samplesOfInterest[["humanEmbryo"]] <- list()
samplesOfInterest[["humanEmbryoOocyte"]] <- list()
samplesOfInterest[["humanXueOocyteEarlyEmbryo"]] <- list()

##### read RPKMs
rpkmFiles <- list.files(RPKM_csv_dir,pattern=".RPKMs.csv$", full.names=TRUE)
rpkms <- lapply(rpkmFiles, read.csv) 
names(rpkms) <- gsub(paste(RPKM_csv_dir,"/",sep=""),"",rpkmFiles)

temp <- as.data.frame(t(sapply(rpkms, dim)))
rownames(temp) <- sapply(strsplit(rownames(temp), "\\."), "[[", 7)
colnames(temp) <- c("numGenes","numColumns")
temp[,"numSamples"] <- (temp[,"numColumns"]-13)/2
temp
# looks good
rm(temp)

head(rpkms[[1]], 1)
# looks good

rpkms_metadata <- data.frame(dataset=names(rpkms))
rpkms_metadata[,"shortName"] <- sapply(strsplit(rpkms_metadata[,"dataset"], "\\."), "[[", 7)
rpkms_metadata[,"species"] <- sapply(strsplit(rpkms_metadata[,"shortName"], "_"), "[[", 1)
rpkms_metadata[,"numSamples"] <- sapply(rpkms, function(x) {  length(grep("_RPKM",colnames(x))) })

### plots, 1 page for each dataset
for (i in 1:dim(rpkms_metadata)[1]) {
  thisDataset <- rpkms_metadata[i,"dataset"]
  thisDatasetShortName <- rpkms_metadata[i,"shortName"]
  thisPlotFile <- paste(plotsDirRaw,"/",gsub("csv$", "pdf", thisDataset),sep="")
  if(file.exists(thisPlotFile)) {next}
  cat("## plotting dataset", thisDatasetShortName, "\n")
  plotRPKMsEachGene(rpkms[[thisDataset]], 
                    plotFile=thisPlotFile, 
                    pageTitles=paste(thisDatasetShortName,"\n\n", thisDataset, "\n", workingDir),
                    textToStripFromSampleNames=paste(rpkms_metadata[i,"species"],"_",sep="") )
  rm(thisDataset,thisDatasetShortName,thisPlotFile)
}
rm(i)




####### special datasets that I want to plot a different way

###### chicken_7organs_development

# make a ggplot-friendly dataframe
tempSamples <- grep("_RPKM",colnames(rpkms[["Chick_genes.blat_galGal6.besthits.psl.bed.counts.chicken_7organs_development.split.multicov.RPKMs.csv"]]), value=TRUE)
chicken7_gg <- data.frame( sample=tempSamples, 
                           rpkm=as.numeric(rpkms[["Chick_genes.blat_galGal6.besthits.psl.bed.counts.chicken_7organs_development.split.multicov.RPKMs.csv"]][,tempSamples]))
rm(tempSamples)

chicken7_gg[,"tissue"] <- sapply(strsplit(chicken7_gg[,"sample"], "_"), "[[", 1)
chicken7_gg[,"devStage"] <- sapply(strsplit(chicken7_gg[,"sample"], "_"), "[[", 2)
chicken7_gg[,"devStage"] <- factor(chicken7_gg[,"devStage"], 
                                   levels=c("e10","e12","e14","e17", "p0","p7","p35","p70","p155"))
chicken7_gg[,"sex"] <- sapply(strsplit(chicken7_gg[,"sample"], "_"), "[[", 3)

table(chicken7_gg[,"tissue"])
# Brain Cerebellum      Heart     Kidney      Liver      Ovary     Testis 
#    36         36         37         36         36         18         18 

table(chicken7_gg[,"devStage"])
# e10  e12  e14  e17   p0   p7  p35  p70 p155 
#  25   24   24   24   24   24   24   24   24 

table(chicken7_gg[,"sex"])
#   F   M 
# 108 109 

pdf(width=11, height=15,
  file=paste(plotsDirRaw,"/Chick_genes.blat_galGal6.besthits.psl.bed.counts.chicken_7organs_development.split.multicov.RPKMs.pdf",sep=""))
ggplot(chicken7_gg, aes(x=devStage, y=rpkm)) +
  facet_grid(rows=vars(tissue),cols=vars(sex)) +
  geom_col() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

##### horse FAANG


# make a ggplot-friendly dataframe
tempSamples <- grep("_RPKM",colnames(rpkms[["Hors_genes.blat_equCab3.besthits.psl.bed.counts.horse_FAANG_tissueSurvey.split.multicov.RPKMs.csv"]]), value=TRUE)
tempGenes <-rpkms[["Hors_genes.blat_equCab3.besthits.psl.bed.counts.horse_FAANG_tissueSurvey.split.multicov.RPKMs.csv"]][,"name"]


horseFAANG_gg <- data.frame( rpkm=t(rpkms[["Hors_genes.blat_equCab3.besthits.psl.bed.counts.horse_FAANG_tissueSurvey.split.multicov.RPKMs.csv"]][,tempSamples]))
colnames(horseFAANG_gg) <- rpkms[["Hors_genes.blat_equCab3.besthits.psl.bed.counts.horse_FAANG_tissueSurvey.split.multicov.RPKMs.csv"]][,"name"]  

## make this a long table not a wide table
horseFAANG_gg_long <- pivot_longer(horseFAANG_gg, cols=everything(), names_to = "Gene", values_to = "RPKM")
horseFAANG_gg_long[,"sample"] <- rep(rownames(horseFAANG_gg), each=dim(horseFAANG_gg)[2])
horseFAANG_gg_long[,"source"] <- sapply(strsplit(pull(horseFAANG_gg_long, sample) , "_"), "[[", 1)
table(horseFAANG_gg_long[,"source"] )
# cells1 cells2  mare1  mare2 
#    40     40    480    480 





horseFAANG_gg_long[,"tissue"] <- gsub("_RPKM","",pull(horseFAANG_gg_long,sample))

horseFAANG_gg_long[,"tissue"] <- gsub("^cells[12]_","",pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("^mare[12]_","",pull(horseFAANG_gg_long, tissue))
## special treatment for some tissues so they are more obviously grouped and easier to understand
horseFAANG_gg_long[,"tissue"] <- gsub("left_ovary","ovary_left", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("right_ovary","ovary_right", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("cerebellar","brain_cerebellar", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("frontal_cortex","brain_frontal_cortex", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("occipital_cortex","brain_occipital_cortex", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("parietal_cortex","brain_parietal_cortex", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("hypothalamus","brain_hypothalamus", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("retina","eye_retina", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("cornea","eye_cornea", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("cardiac_muscle_of_left_atrium","heart_left_atrium", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("mitral_valve","heart_mitral_valve", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("longissimus_thoracis_muscle","muscle_longissimus_thoracis", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("ventral_lateral_sacrocaudal_muscle","muscle_ventral_lateral_sacrocaudal", pull(horseFAANG_gg_long, tissue))
horseFAANG_gg_long[,"tissue"] <- gsub("gluteal_muscle","muscle_gluteal", pull(horseFAANG_gg_long, tissue))

horseFAANG_gg_long[,"tissueNoRep"] <- gsub("_rep[12]","", pull(horseFAANG_gg_long, tissue))
table(horseFAANG_gg_long[,"tissueNoRep"] )


pdf(width=11, height=15,
    file=paste(plotsDirRaw,"/Hors_genes.blat_equCab3.besthits.psl.bed.counts.horse_FAANG_tissueSurvey.split.multicov.RPKMs.ggplot.pdf",sep=""))
ggplot(horseFAANG_gg_long, aes_string(x="tissueNoRep", y="RPKM")) +
  facet_grid(rows=vars(Gene)) +
  geom_jitter(size=0.25, width=0.1, height=0) + 
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()



###### human oogenesis - BEFORE I trimmed the reads and remapped

human_oogenesis <- rpkms[["Hum_genes.blat_hg38.besthits.psl.bed.counts.human_oogenesis.split.multicov.RPKMs.csv"]]
tempCols <- grep("RPKM", colnames(human_oogenesis))

human_oogenesis_gg <- t(human_oogenesis[,tempCols])
colnames(human_oogenesis_gg) <- human_oogenesis[,"name"]

## make this a long table not a wide table
human_oogenesis_gg_long <- pivot_longer(as.data.frame(human_oogenesis_gg), cols=everything(), names_to = "Gene", values_to = "RPKM")
human_oogenesis_gg_long[,"sample"] <- rep(rownames(human_oogenesis_gg), each=dim(human_oogenesis_gg)[2])
human_oogenesis_gg_long[,"sample"] <- gsub("Preovu","preovu",pull(human_oogenesis_gg_long, sample))
human_oogenesis_gg_long[,"sampleType"] <- gsub("_RPKM","",pull(human_oogenesis_gg_long, sample))
human_oogenesis_gg_long[,"sampleType"] <- gsub("_\\d+","",pull(human_oogenesis_gg_long, sampleType))
human_oogenesis_gg_long[,"sampleType"] <- factor(pull(human_oogenesis_gg_long, sampleType), 
      levels=c("GC_primordial", "GC_primary", "GC_secondary", "GC_antral", "GC_preovulatory",
               "Oocyte_primordial", "Oocyte_primary", "Oocyte_secondary", "Oocyte_antral"))
human_oogenesis_gg_long[,"cellType"] <- sapply(strsplit(as.character(pull(human_oogenesis_gg_long, sampleType)), "_"), "[[", 1)
human_oogenesis_gg_long[,"stage"] <- sapply(strsplit(as.character(pull(human_oogenesis_gg_long, sampleType)), "_"), "[[", 2)
human_oogenesis_gg_long[,"stage"] <- factor(pull(human_oogenesis_gg_long, stage), 
                  levels=c("primordial", "primary", "secondary", "antral", "preovulatory"))

pdf(width=11, height=15,
    file=paste(plotsDirRaw,"/Hum_genes.blat_hg38.besthits.psl.bed.counts.human_oogenesis.split.multicov.RPKMs.ggplot.pdf",sep=""))
ggplot(human_oogenesis_gg_long, aes_string(x="sampleType", y="RPKM")) +
  facet_grid(rows=vars(Gene), scales="free_y") +
  geom_violin(scale = "width") + 
  stat_summary(fun=median, geom="point", size=2, color="red") +
  geom_jitter(size=0.25, width=0.1, height=0) + 
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

## try plots that look more like the bar graphs I used

human_oogenesis_gg_long_df <- as.data.frame(human_oogenesis_gg_long)

human_oogenesis_gg_long_df[,"geneShort"] <- gsub("_Hum", "", human_oogenesis_gg_long_df[,"Gene"])
human_oogenesis_gg_long_df[,"geneShort"] <- gsub("H2Bopp8C", "H2B.N", human_oogenesis_gg_long_df[,"geneShort"])
human_oogenesis_gg_long_df[,"geneShort"] <- gsub("-pseudo", "", human_oogenesis_gg_long_df[,"geneShort"])



human_oogenesis_gg_long_df[,"geneShort"] <- factor(human_oogenesis_gg_long_df[,"geneShort"], 
          levels= intersect(preferredGeneOrder, 
                      unique(human_oogenesis_gg_long_df[,"geneShort"]))   )
human_oogenesis_gg_long_summary <- human_oogenesis_gg_long_df %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

## plots like I had for the other datasets
human_oogenesis_plot <- ggplot(human_oogenesis_gg_long_summary, 
                               aes(x=sampleType, y=my_median, fill=geneShort, 
                                   ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  #guides(fill=guide_legend(title="gene")) +
  genesFillScale +
  labs(x="", y="RPKM", title="Human oogenesis") +
  theme_classic() +   
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))

human_oogenesis_plot

## violin plots 
human_oogenesis_violinplot <- ggplot(human_oogenesis_gg_long_df, 
                               aes(x=sampleType, y=RPKM, fill=geneShort, color=geneShort)) +
  geom_violin(scale = "width") + 
  stat_summary(fun=median, geom="point", size=0.5, 
               color="red", position=position_dodge(0.9)) +
  genesFillScale + genesColorScale +
  guides(fill=guide_legend(override.aes = list(shape = NA))) +
  labs(x="", y="RPKM", title="Human oogenesis") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
human_oogenesis_violinplot



########## put all those plots on one page for output:
pdf(file=paste(plotsDirPretty,"/human_oogenesisPlots_2021_Apr28.pdf",sep=""), 
    width=11, height=5)
ggarrange(human_oogenesis_violinplot, 
          human_oogenesis_violinplot+coord_cartesian(ylim=c(0,100)), 
          ncol=2, nrow=1)
dev.off()



###### human oogenesis - AFTER I trimmed the reads and remapped

human_oogenesisTrimReads <- rpkms[["Hum_genes.blat_hg38.besthits.psl.bed.counts.human_oogenesisTrimReads.split.multicov.RPKMs.csv"]]
tempCols <- grep("RPKM", colnames(human_oogenesisTrimReads))

human_oogenesisTrimReads_gg <- t(human_oogenesisTrimReads[,tempCols])
colnames(human_oogenesisTrimReads_gg) <- human_oogenesisTrimReads[,"name"]

## make this a long table not a wide table
human_oogenesisTrimReads_gg_long <- pivot_longer(as.data.frame(human_oogenesisTrimReads_gg), cols=everything(), names_to = "Gene", values_to = "RPKM")
human_oogenesisTrimReads_gg_long[,"sample"] <- rep(rownames(human_oogenesisTrimReads_gg), each=dim(human_oogenesisTrimReads_gg)[2])
human_oogenesisTrimReads_gg_long[,"sample"] <- gsub("Preovu","preovu",pull(human_oogenesisTrimReads_gg_long, sample))
human_oogenesisTrimReads_gg_long[,"sampleType"] <- gsub("_RPKM","",pull(human_oogenesisTrimReads_gg_long, sample))
human_oogenesisTrimReads_gg_long[,"sampleType"] <- gsub("_\\d+","",pull(human_oogenesisTrimReads_gg_long, sampleType))
human_oogenesisTrimReads_gg_long[,"sampleType"] <- factor(pull(human_oogenesisTrimReads_gg_long, sampleType), 
  levels=c("GC_primordial", "GC_primary", "GC_secondary", "GC_antral", "GC_preovulatory",
            "Oocyte_primordial", "Oocyte_primary", "Oocyte_secondary", "Oocyte_antral"))
human_oogenesisTrimReads_gg_long[,"cellType"] <- sapply(strsplit(as.character(pull(human_oogenesisTrimReads_gg_long, sampleType)), "_"), "[[", 1)
human_oogenesisTrimReads_gg_long[,"stage"] <- sapply(strsplit(as.character(pull(human_oogenesisTrimReads_gg_long, sampleType)), "_"), "[[", 2)
human_oogenesisTrimReads_gg_long[,"stage"] <- factor(pull(human_oogenesisTrimReads_gg_long, stage), 
     levels=c("primordial", "primary", "secondary", "antral", "preovulatory"))





### violin plots, kinda ugly, not using
pdf(width=11, height=15,
    file=paste(plotsDirRaw,"/Hum_genes.blat_hg38.besthits.psl.bed.counts.human_oogenesisTrimReads.split.multicov.RPKMs.ggplot.pdf",sep=""))
ggplot(human_oogenesisTrimReads_gg_long, aes_string(x="sampleType", y="RPKM")) +
  facet_grid(rows=vars(Gene), scales="free_y") +
  geom_violin(scale = "width") + 
  stat_summary(fun=median, geom="point", size=2, color="red") +
  geom_jitter(size=0.25, width=0.1, height=0) + 
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

## try plots that look more like the bar graphs I used
human_oogenesisTrimReads_gg_long_df <- as.data.frame(human_oogenesisTrimReads_gg_long)

human_oogenesisTrimReads_gg_long_df[,"geneShort"] <- gsub("_Hum", "", human_oogenesisTrimReads_gg_long_df[,"Gene"])
human_oogenesisTrimReads_gg_long_df[,"geneShort"] <- gsub("H2Bopp8C", "H2B.N", human_oogenesisTrimReads_gg_long_df[,"geneShort"])
human_oogenesisTrimReads_gg_long_df[,"geneShort"] <- gsub("-pseudo", "", human_oogenesisTrimReads_gg_long_df[,"geneShort"])

human_oogenesisTrimReads_gg_long_df[,"geneShort"] <- factor(human_oogenesisTrimReads_gg_long_df[,"geneShort"], 
         levels= intersect(preferredGeneOrder, 
                           unique(human_oogenesisTrimReads_gg_long_df[,"geneShort"]))   )

human_oogenesisTrimReads_gg_long_summary <- human_oogenesisTrimReads_gg_long_df %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

## plots like I had for the other datasets
human_oogenesisTrimReads_plot <- ggplot(human_oogenesisTrimReads_gg_long_summary, 
                               aes(x=sampleType, y=my_median, fill=geneShort, 
                                   ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Human oogenesis, trimmed reads") +
  theme_classic() +   
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))

human_oogenesisTrimReads_plot

## violin plots 
human_oogenesisTrimReads_violinplot <- ggplot(human_oogenesisTrimReads_gg_long_df, 
                                     aes(x=sampleType, y=RPKM, fill=geneShort, color=geneShort)) +
  geom_violin(scale = "width") + 
  stat_summary(fun=median, geom="point", size=0.5, 
               color="red", position=position_dodge(0.9)) +
  genesFillScale + genesColorScale +
  guides(fill=guide_legend(title="gene", override.aes = list(shape = NA)), 
         color=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Human oogenesis, trimmed reads") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
human_oogenesisTrimReads_violinplot



########## put all those plots on one page for output:
pdf(file=paste(plotsDirPretty,"/human_oogenesisTrimReadsPlots_2021_Apr28.pdf",sep=""), 
    width=11, height=5)
ggarrange(human_oogenesisTrimReads_violinplot, 
          human_oogenesisTrimReads_violinplot+coord_cartesian(ylim=c(0,100)), 
          ncol=2, nrow=1)
dev.off()



##########for a figure: plot just ancestral copies of each variant

######## dog  
# use ovary from dog_transcriptome2 not from dog_selectedTissues)
# we do not have embryo


samplesOfInterest[["dog"]] [["dog_selectedTissues"]] <- c(
  "brain_rep1", "brain_rep2", "brain_rep3", 
  "liver_rep1", "liver_rep2", "liver_rep3", 
  "kidney_rep1", "kidney_rep2", "kidney_rep3", 
  "heart_rep1", "heart_rep2", "heart_rep3", 
  "testis_rep1", "testis_rep2", "testis_rep3")
samplesOfInterest[["dog"]] [["dog_transcriptome2"]] <- c("ovary")

# get RPKMs for just the tissues of interest:
dogSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="dog")
# select only some genes
dogSelectedTissues <- dogSelectedTissues[genesOfInterest[["dog"]],]

# convert to long format
dogSelectedTissues_long <- convertToLongFormat(dogSelectedTissues)
dogSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d","",
             gsub("_RPKM","",dogSelectedTissues_long[,"sample"]))
dogSelectedTissues_long[,"sampleType"] <- factor(
  dogSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(dogSelectedTissues_long[,"sampleType"]))) 

dogSelectedTissues_long[,"geneShort"] <- factor(sapply(strsplit(
  dogSelectedTissues_long[,"Gene"], "_Dog"),"[[", 1), levels=preferredGeneOrder)

dogSelectedTissues_long_summary <- dogSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data - this looks good:
dog_plot <- ggplot(dogSelectedTissues_long_summary, 
       aes(x=sampleType, y=my_median, fill=geneShort, 
           ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Dog") +
  theme_classic()
dog_plot



######## pig  


samplesOfInterest[["pig"]] [["pig_selectedTissues"]] <- c("brainFrontalLobe", "liver", "kidney", "heart", "testis_Duroc", "testis_LargeWhite", "ovary")
# I am deliberately not using the 3rd testis sample - not positive for anything (?) (xx would be good to check a control gene)

# get RPKMs for just the tissues of interest:
pigSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="pig")

# select only some genes
pigSelectedTissues <- pigSelectedTissues[genesOfInterest[["pig"]],]

# convert to long format
pigSelectedTissues_long <- convertToLongFormat(pigSelectedTissues)
pigSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",pigSelectedTissues_long[,"sample"])
pigSelectedTissues_long[,"sampleType"] <- gsub("brainFrontalLobe","brain",pigSelectedTissues_long[,"sampleType"])
pigSelectedTissues_long[,"sampleType"] <- gsub("testis_Duroc","testis",pigSelectedTissues_long[,"sampleType"])
pigSelectedTissues_long[,"sampleType"] <- gsub("testis_LargeWhite","testis",pigSelectedTissues_long[,"sampleType"])


pigSelectedTissues_long[,"sampleType"] <- factor(
  pigSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(pigSelectedTissues_long[,"sampleType"]))) 

pigSelectedTissues_long[,"geneShort"] <- sapply(strsplit(
  pigSelectedTissues_long[,"Gene"], "_Pig"),"[[", 1)

pigSelectedTissues_long[,"geneShort"] <- factor(pigSelectedTissues_long[,"geneShort"], 
                                                levels= intersect(preferredGeneOrder, 
                    unique(pigSelectedTissues_long[,"geneShort"]))   )

pigSelectedTissues_long_summary <- pigSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
pig_plot <- ggplot(pigSelectedTissues_long_summary, 
                   aes(x=sampleType, y=my_median, fill=geneShort, 
                       ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Pig") +
  theme_classic()
pig_plot





######## mouse  


samplesOfInterest[["mouse"]] [["mouse_Brawand"]] <- c(
  "mouse_Brain1", "mouse_Brain2", "mouse_Brain3", 
  "mouse_Brain4", "mouse_Brain5", "mouse_Brain6", 
  "mouse_Liver1", "mouse_Liver2", "mouse_Liver3", 
  "mouse_Kidney1", "mouse_Kidney2", "mouse_Kidney3", 
  "mouse_Heart1", "mouse_Heart2", "mouse_Heart3", 
  "mouse_Testis1", "mouse_Testis2")
samplesOfInterest[["mouse"]] [["mouse_oocyte"]] <- c("oocytes_nonGrowing", "oocytes_Growing_rep1", "oocytes_Growing_rep2", "oocytes_GV_rep1", "oocytes_GV_rep2", "oocytes_GV_rep3", "oocytes_GV_rep4", "oocytes_metII_rep1", "oocytes_metII_rep2", "oocytes_metII_rep3", "oocytes_E18.5")

# get RPKMs for just the tissues of interest:
mouseSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="mouse")

# select only some genes
mouseSelectedTissues <- mouseSelectedTissues[genesOfInterest[["mouse"]],]

# convert to long format
mouseSelectedTissues_long <- convertToLongFormat(mouseSelectedTissues)
mouseSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",mouseSelectedTissues_long[,"sample"])
mouseSelectedTissues_long[,"sampleType"] <- gsub("mouse_","",mouseSelectedTissues_long[,"sampleType"])
mouseSelectedTissues_long[,"sampleType"] <- gsub("\\d","",mouseSelectedTissues_long[,"sampleType"])
mouseSelectedTissues_long[,"sampleType"] <- sapply(strsplit(mouseSelectedTissues_long[,"sampleType"], "_"), "[[", 1)
mouseSelectedTissues_long[,"sampleType"] <- tolower(mouseSelectedTissues_long[,"sampleType"] )
mouseSelectedTissues_long[,"sampleType"] <- factor(
  mouseSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(mouseSelectedTissues_long[,"sampleType"]))) 

mouseSelectedTissues_long[,"geneShort"] <- sapply(strsplit(
  mouseSelectedTissues_long[,"Gene"], "_Mouse"),"[[", 1)
mouseSelectedTissues_long[,"geneShort"] <- factor(mouseSelectedTissues_long[,"geneShort"], 
                                                levels= intersect(preferredGeneOrder, 
                                                                  unique(mouseSelectedTissues_long[,"geneShort"]))   )

mouseSelectedTissues_long_summary <- mouseSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
mouse_plot <- ggplot(mouseSelectedTissues_long_summary, 
                   aes(x=sampleType, y=my_median, fill=geneShort, 
                       ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Mouse") +
  theme_classic()
mouse_plot


##########  mouse_PGCs

samplesOfInterest[["mouse"]] [["mouse_PGCs"]] <- gsub("_reads","",grep("reads", colnames(rpkms[["Mouse_genes.blat_mouse_Dec2011.besthits.psl.bed.counts.mouse_PGCs.split.multicov.RPKMs.csv"]]), value = TRUE))
  
dat <- rpkms[["Mouse_genes.blat_mouse_Dec2011.besthits.psl.bed.counts.mouse_PGCs.split.multicov.RPKMs.csv"]] %>% 
  select(name, ends_with("RPKM"))

dat_long <- dat %>% 
  pivot_longer(cols=-name, 
               names_sep="_",
               names_to=c("devStage","sex",NA,NA),
               values_to="RPKM")  %>% # NA means discard that bit
  mutate(geneName=gsub("_Mouse","",name)) %>% 
  mutate(geneName=factor(geneName, levels=intersect(preferredGeneOrder, geneName)) ) 


### plot
mouse_PGC_plot <- dat_long %>% 
  ggplot(aes(x=devStage, y=RPKM)) +
  geom_point(aes(fill=geneName, colour=geneName)) +
  stat_summary(colour="dark gray", #aes(colour=geneName), 
               fun = median, geom = "crossbar", 
               width = 0.5, fatten=1.5, show.legend=FALSE) +
  facet_grid(rows=vars(sex), cols=vars(geneName)) +
  labs(x="", y="RPKM", title="Mouse PGCs") +
  genesFillScale + 
  genesColorScale +
  guides(colour=guide_legend(title="gene"),fill=guide_legend(title="gene")) +
  #theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(mouse_PGC_plot, device="pdf", file="aa_plots/mouse_PGCPlots_2021_Apr19.pdf")

######## human  


samplesOfInterest[["human"]] [["human_Brawand"]] <- c("human_Brain_frontal_cortex1", "human_Brain_frontal_cortex2", "human_Brain_prefrontal_cortex", "human_Brain_prefrontal_cortex_Male1", "human_Brain_prefrontal_cortex_Male2", "human_Brain_temporal_lobe_Male", "human_Liver1", "human_Liver2", "human_Liver3", "human_Kidney_Female", "human_Kidney_Male1", "human_Kidney_Male2", "human_Heart_Female", "human_Heart_Male1", "human_Heart_Male2", "human_Heart_Male3", "human_Testis1", "human_Testis2")
samplesOfInterest[["human"]] [["human_ovary"]] <- c("ovary_fetal_20wk", "ovary_fetal_19wk", "ovary_fetal_16wk")




# get RPKMs for just the tissues of interest:
humanSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="human")

# select only some genes
humanSelectedTissues <- humanSelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanSelectedTissues_long <- convertToLongFormat(humanSelectedTissues)
humanSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanSelectedTissues_long[,"sample"])
humanSelectedTissues_long[,"sampleType"] <- gsub("human_","",humanSelectedTissues_long[,"sampleType"])
humanSelectedTissues_long[,"sampleType"] <- gsub("\\d","",humanSelectedTissues_long[,"sampleType"])
humanSelectedTissues_long[,"sampleType"] <- sapply(strsplit(humanSelectedTissues_long[,"sampleType"], "_"), "[[", 1)
humanSelectedTissues_long[,"sampleType"] <- tolower(humanSelectedTissues_long[,"sampleType"] )
table(humanSelectedTissues_long[,"sampleType"] )

humanSelectedTissues_long[,"sampleType"] <- factor(
  humanSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(humanSelectedTissues_long[,"sampleType"]))) 

humanSelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanSelectedTissues_long[,"Gene"])
humanSelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanSelectedTissues_long[,"geneShort"])

humanSelectedTissues_long[,"geneShort"] <- factor(humanSelectedTissues_long[,"geneShort"], 
                                                  levels= intersect(preferredGeneOrder, 
                                                                    unique(humanSelectedTissues_long[,"geneShort"]))   )

humanSelectedTissues_long_summary <- humanSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
human_plot <- ggplot(humanSelectedTissues_long_summary, 
                     aes(x=sampleType, y=my_median, fill=geneShort, 
                         ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Human") +
  theme_classic()
human_plot




######## chicken  


samplesOfInterest[["chicken"]] [["chicken_21tissues"]] <- c("Cerebellum", "Liver", "Kidney", "Heart_Muscle", "Ovary")
samplesOfInterest[["chicken"]] [["chicken_Brawand"]] <- c("chicken_Testis_Male1", "chicken_Testis_Male2", "chicken_Testis_Male3")


# get RPKMs for just the tissues of interest:
chickenSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="chicken")

# select only some genes
chickenSelectedTissues <- chickenSelectedTissues[genesOfInterest[["chicken"]],]

# convert to long format
chickenSelectedTissues_long <- convertToLongFormat(chickenSelectedTissues)
chickenSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",chickenSelectedTissues_long[,"sample"])
chickenSelectedTissues_long[,"sampleType"] <- gsub("Cerebellum","brain",chickenSelectedTissues_long[,"sampleType"])
chickenSelectedTissues_long[,"sampleType"] <- gsub("chicken_","",chickenSelectedTissues_long[,"sampleType"])
chickenSelectedTissues_long[,"sampleType"] <- sapply(strsplit(chickenSelectedTissues_long[,"sampleType"], "_"), "[[", 1)
chickenSelectedTissues_long[,"sampleType"] <- tolower(chickenSelectedTissues_long[,"sampleType"] )
table(chickenSelectedTissues_long[,"sampleType"])

chickenSelectedTissues_long[,"sampleType"] <- factor(
  chickenSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(chickenSelectedTissues_long[,"sampleType"]))) 

chickenSelectedTissues_long[,"geneShort"] <- sapply(strsplit(
  chickenSelectedTissues_long[,"Gene"], "_Chick"),"[[", 1)

chickenSelectedTissues_long[,"geneShort"] <- factor(chickenSelectedTissues_long[,"geneShort"], 
                                                  levels= intersect(preferredGeneOrder, 
                                                                    unique(chickenSelectedTissues_long[,"geneShort"]))   )

chickenSelectedTissues_long_summary <- chickenSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data: (no error bars, only one sample per tissue)
chicken_plot <- ggplot(chickenSelectedTissues_long_summary, 
                     aes(x=sampleType, y=my_median, fill=geneShort, 
                         ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  #geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="chicken") +
  theme_classic()
chicken_plot




######## opossum  


samplesOfInterest[["opossum"]] [["opossum_Brawand"]] <- c("opossum_Brain_Female1", "opossum_Brain_Female2", "opossum_Brain_Male", "opossum_Liver_Female", "opossum_Liver_Male", "opossum_Kidney_Female", "opossum_Kidney_Male", "opossum_Heart_Female1", "opossum_Heart_Female2", "opossum_Heart_Male1", "opossum_Heart_Male2", "opossum_Testis_Male1", "opossum_Testis_Male2")
samplesOfInterest[["opossum"]] [["opossum_ovary"]] <- c("ovary1", "ovary2", "ovary3", "ovary4", "ovary5", "ovary6", "ovary7")



# get RPKMs for just the tissues of interest:
opossumSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="opossum")

# select only some genes
opossumSelectedTissues <- opossumSelectedTissues[genesOfInterest[["opossum"]],]

# convert to long format
opossumSelectedTissues_long <- convertToLongFormat(opossumSelectedTissues)
opossumSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",opossumSelectedTissues_long[,"sample"])
opossumSelectedTissues_long[,"sampleType"] <- gsub("opossum_","",opossumSelectedTissues_long[,"sampleType"])
opossumSelectedTissues_long[,"sampleType"] <- gsub("\\d","",opossumSelectedTissues_long[,"sampleType"])
opossumSelectedTissues_long[,"sampleType"] <- sapply(strsplit(opossumSelectedTissues_long[,"sampleType"], "_"), "[[", 1)
opossumSelectedTissues_long[,"sampleType"] <- tolower(opossumSelectedTissues_long[,"sampleType"] )

opossumSelectedTissues_long[,"sampleType"] <- factor(
  opossumSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(opossumSelectedTissues_long[,"sampleType"]))) 

opossumSelectedTissues_long[,"geneShort"] <- sapply(strsplit(
  opossumSelectedTissues_long[,"Gene"], "_Oppos"),"[[", 1)

opossumSelectedTissues_long[,"geneShort"] <- factor(opossumSelectedTissues_long[,"geneShort"], 
                                                    levels= intersect(preferredGeneOrder, 
                                                                      unique(opossumSelectedTissues_long[,"geneShort"]))   )

opossumSelectedTissues_long_summary <- opossumSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
opossum_plot <- ggplot(opossumSelectedTissues_long_summary, 
                       aes(x=sampleType, y=my_median, fill=geneShort, 
                           ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Opossum") +
  theme_classic()
opossum_plot




######## macaque  


samplesOfInterest[["macaque"]] [["macaque_Brawand"]] <- c("macaque_Brain_Female", "macaque_Brain_Male", "macaque_Brain_prefrontal_cortex_Male", "macaque_Liver_Female", "macaque_Liver_Male1", "macaque_Liver_Male2", "macaque_Kidney_Female", "macaque_Kidney_Male", "macaque_Heart_Female", "macaque_Heart_Male", "macaque_Testis_Male1", "macaque_Testis_Male2")



# get RPKMs for just the tissues of interest:
macaqueSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="macaque")

# select only some genes
macaqueSelectedTissues <- macaqueSelectedTissues[genesOfInterest[["macaque"]],]

# convert to long format
macaqueSelectedTissues_long <- convertToLongFormat(macaqueSelectedTissues)
macaqueSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",macaqueSelectedTissues_long[,"sample"])
macaqueSelectedTissues_long[,"sampleType"] <- gsub("macaque_","",macaqueSelectedTissues_long[,"sampleType"])
macaqueSelectedTissues_long[,"sampleType"] <- sapply(strsplit(macaqueSelectedTissues_long[,"sampleType"], "_"), "[[", 1)
macaqueSelectedTissues_long[,"sampleType"] <- tolower(macaqueSelectedTissues_long[,"sampleType"] )

macaqueSelectedTissues_long[,"sampleType"] <- factor(
  macaqueSelectedTissues_long[,"sampleType"],
  levels=intersect(preferredSampleTypeOrder, 
                   unique(macaqueSelectedTissues_long[,"sampleType"]))) 

macaqueSelectedTissues_long[,"geneShort"] <- sapply(strsplit(
  macaqueSelectedTissues_long[,"Gene"], "_Macaque"),"[[", 1)

macaqueSelectedTissues_long[,"geneShort"] <- factor(macaqueSelectedTissues_long[,"geneShort"], 
                                                    levels= intersect(preferredGeneOrder, 
                                                                      unique(macaqueSelectedTissues_long[,"geneShort"]))   )

macaqueSelectedTissues_long_summary <- macaqueSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
macaque_plot <- ggplot(macaqueSelectedTissues_long_summary, 
                       aes(x=sampleType, y=my_median, fill=geneShort, 
                           ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="Macaque") +
  theme_classic()
macaque_plot




########## put all those plots on one page for output:
pdf(file=paste(plotsDirPretty,"/combinedPlots_2021_Mar24.pdf",sep=""), width=11, height=25)
ggarrange(human_plot, human_plot+coord_cartesian(ylim=c(0,4)), 
          macaque_plot, macaque_plot+coord_cartesian(ylim=c(0,2)),
          mouse_plot, mouse_plot+coord_cartesian(ylim=c(0,30)), 
          pig_plot, pig_plot+coord_cartesian(ylim=c(0,30)), 
          dog_plot, dog_plot+coord_cartesian(ylim=c(0,30)),
          opossum_plot, opossum_plot+coord_cartesian(ylim=c(0,100)), 
          chicken_plot, chicken_plot+coord_cartesian(ylim=c(0,1)), 
          ncol=2, nrow=7)
dev.off()

# also show with square-root transformation
pdf(file=paste(plotsDirPretty,"/combinedPlots_2021_Mar24_sqrtTrans.pdf",sep=""), width=7, height=25)
ggarrange(human_plot + scale_y_sqrt(), 
          macaque_plot + scale_y_sqrt(), 
          mouse_plot + scale_y_sqrt(), 
          pig_plot + scale_y_sqrt(), 
          dog_plot + scale_y_sqrt(), 
          opossum_plot + scale_y_sqrt(), 
          chicken_plot + scale_y_sqrt(), 
          ncol=1, nrow=7)
dev.off()






######## human spermatogenesis 

# I took the spermatogenesis samples from the human_testis_sperm datset and split them, because the Whitehead data ("study 1") gives me RPKMs on a VERY different scale than the other data, from a paper by Jan et al (2017) ("study 2"), see further down

### study 1

samplesOfInterest[["humanSpermatogenesisStudy1"]] [["human_testis_sperm"]] <- 
  c("human1_pachytene_sperm", "human1_roundspermatid_sperm", "human2_pachytene_sperm",
    "human2_roundspermatid_sperm", "human3_pachytene_sperm", "human3_roundspermatid_sperm")

# get RPKMs for just the tissues of interest:
humanSpermatogenesisStudy1SelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="humanSpermatogenesisStudy1")

# select only some genes
humanSpermatogenesisStudy1SelectedTissues <- humanSpermatogenesisStudy1SelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanSpermatogenesisStudy1SelectedTissues_long <- convertToLongFormat(humanSpermatogenesisStudy1SelectedTissues)
humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanSpermatogenesisStudy1SelectedTissues_long[,"sample"])
humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"] <- gsub("human\\d_","",humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"])
humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"] <- gsub("roundspermatid_sperm","round_spermatids",humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"])
humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"] <- gsub("^pachytene_sperm","pachytene_spermatocytes",humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"])

# order the sperm types (see fig 1 in Jan et al paper in /fh/fast/malik_h/grp/public_databases/NCBI/SRA/data/mammalian_expression_profiles/human/human_spermatogenesis)
humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"] <- factor(
  humanSpermatogenesisStudy1SelectedTissues_long[,"sampleType"],
  levels=c("pachytene_spermatocytes", 
           "round_spermatids")) 

humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanSpermatogenesisStudy1SelectedTissues_long[,"Gene"])
humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"])

humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"] <- factor(
  humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanSpermatogenesisStudy1SelectedTissues_long[,"geneShort"]))   )

humanSpermatogenesisStudy1SelectedTissues_long_summary <- humanSpermatogenesisStudy1SelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
humanSpermatogenesisStudy1_plot <- ggplot(humanSpermatogenesisStudy1SelectedTissues_long_summary, 
                       aes(x=sampleType, y=my_median, fill=geneShort, 
                           ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human spermatogenesis (Whitehead)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))

#### human spermatogenesis study 2
# I initially got 1 replicate of each sample type
# updating to show all 34 samples

samplesOfInterest[["humanSpermatogenesisStudy2"]] [["human_spermatogenesis"]] <- 
  c("dark_spermatogonia_1", "dark_spermatogonia_2", "dark_spermatogonia_3", 
    "dark_spermatogonia_4", "dark_spermatogonia_5", "dark_spermatogonia_6", 
    "early_pachytene_spermatocytes_1", "early_pachytene_spermatocytes_2", 
    "early_pachytene_spermatocytes_3", "early_pachytene_spermatocytes_4", 
    "early_pachytene_spermatocytes_5", "early_pachytene_spermatocytes_6", 
    "late_pachytene_spermatocytes_1", "late_pachytene_spermatocytes_2", 
    "late_pachytene_spermatocytes_3", "late_pachytene_spermatocytes_4", 
    "late_pachytene_spermatocytes_5", "late_pachytene_spermatocytes_6", 
    "leptotene_zygotene_spermatocytes_1", "leptotene_zygotene_spermatocytes_2", 
    "leptotene_zygotene_spermatocytes_3", "leptotene_zygotene_spermatocytes_4", 
    "leptotene_zygotene_spermatocytes_5", 
    "pale_spermatogonia_1", "pale_spermatogonia_2", "pale_spermatogonia_3", 
    "pale_spermatogonia_4", "pale_spermatogonia_5", 
    "round_spermatids_1", "round_spermatids_2", "round_spermatids_3", 
    "round_spermatids_4", "round_spermatids_5", "round_spermatids_6")

# get RPKMs for just the tissues of interest:
humanSpermatogenesisStudy2SelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, species="humanSpermatogenesisStudy2")

# select only some genes
humanSpermatogenesisStudy2SelectedTissues <- humanSpermatogenesisStudy2SelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanSpermatogenesisStudy2SelectedTissues_long <- convertToLongFormat(humanSpermatogenesisStudy2SelectedTissues)
humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanSpermatogenesisStudy2SelectedTissues_long[,"sample"])
humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"] <- gsub("_\\d$","",humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"])

table(humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"])

# order the sperm types (see fig 1 in Jan et al paper in /fh/fast/malik_h/grp/public_databases/NCBI/SRA/data/mammalian_expression_profiles/human/human_spermatogenesis)
humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"] <- factor(
  humanSpermatogenesisStudy2SelectedTissues_long[,"sampleType"],
  levels=c("dark_spermatogonia", 
           "pale_spermatogonia", 
           "leptotene_zygotene_spermatocytes", 
           "early_pachytene_spermatocytes", 
           "late_pachytene_spermatocytes",
           "round_spermatids")) 


humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanSpermatogenesisStudy2SelectedTissues_long[,"Gene"])
humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"])

humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"] <- factor(
  humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanSpermatogenesisStudy2SelectedTissues_long[,"geneShort"]))   )


humanSpermatogenesisStudy2SelectedTissues_long_summary <- humanSpermatogenesisStudy2SelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data - no error bars for now as there's only 1 sample per type
humanSpermatogenesisStudy2_plot <- ggplot(humanSpermatogenesisStudy2SelectedTissues_long_summary, 
                                          aes(x=sampleType, y=my_median, fill=geneShort, 
                                              ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human spermatogenesis (Jan et al 2017)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))


########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/human_spermatogenesisPlots_2021_Mar25.pdf",sep=""), 
    width=11, height=11)
ggarrange(humanSpermatogenesisStudy2_plot, 
          humanSpermatogenesisStudy2_plot+coord_cartesian(ylim=c(0,10)),
          humanSpermatogenesisStudy1_plot, 
          humanSpermatogenesisStudy1_plot+coord_cartesian(ylim=c(0,10)), 
          ncol=2, nrow=2)
dev.off()





########## put all the 6-tissue plots on one page for output:
pdf(file=paste(plotsDirPretty,"/combinedPlots_2021_Mar23.pdf",sep=""), width=11, height=25)
ggarrange(human_plot, human_plot+coord_cartesian(ylim=c(0,4)), 
          macaque_plot, macaque_plot+coord_cartesian(ylim=c(0,2)),
          mouse_plot, mouse_plot+coord_cartesian(ylim=c(0,30)), 
          pig_plot, pig_plot+coord_cartesian(ylim=c(0,30)), 
          dog_plot, dog_plot+coord_cartesian(ylim=c(0,30)),
          opossum_plot, opossum_plot+coord_cartesian(ylim=c(0,100)), 
          chicken_plot, chicken_plot+coord_cartesian(ylim=c(0,1)), 
          ncol=2, nrow=7)
dev.off()

# also show with square-root transformation
pdf(file=paste(plotsDirPretty,"/combinedPlots_2021_Mar23_sqrtTrans.pdf",sep=""), width=7, height=25)
ggarrange(human_plot + scale_y_sqrt(), 
          macaque_plot + scale_y_sqrt(), 
          mouse_plot + scale_y_sqrt(), 
          pig_plot + scale_y_sqrt(), 
          dog_plot + scale_y_sqrt(), 
          opossum_plot + scale_y_sqrt(), 
          chicken_plot + scale_y_sqrt(), 
          ncol=1, nrow=7)
dev.off()



######## human embryo - 
## xxx I can update this later because I am adding some samples, in human_embryoOocyte dir

samplesOfInterest[["humanEmbryo"]][["human_embryo"]] <- 
  c("embryo_CL_rep1a", "embryo_CL_rep1b", "embryo_CL_rep2a", "embryo_CL_rep2b", 
    "embryo_ICM_rep1a", "embryo_ICM_rep1b", "embryo_ICM_rep2a", "embryo_ICM_rep2b", 
    "embryo_MOR_rep1a", "embryo_MOR_rep1b", "embryo_MOR_rep2a", "embryo_MOR_rep2b", 
    "embryo_PN_rep2a", "embryo_PN_rep2b", 
    "embryo_TROPH_rep1a", "embryo_TROPH_rep1b", "embryo_TROPH_rep2a", "embryo_TROPH_rep2b")

# get RPKMs for just the tissues of interest:
humanEmbryoSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                    species="humanEmbryo")

# select only some genes
humanEmbryoSelectedTissues <- humanEmbryoSelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanEmbryoSelectedTissues_long <- convertToLongFormat(humanEmbryoSelectedTissues)
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoSelectedTissues_long[,"sample"])


# convert to long format
humanEmbryoSelectedTissues_long <- convertToLongFormat(humanEmbryoSelectedTissues)
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoSelectedTissues_long[,"sample"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d[ab]$","",humanEmbryoSelectedTissues_long[,"sampleType"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_MOR","morula", humanEmbryoSelectedTissues_long[,"sampleType"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_CL","2-4-8-cell_embryo",humanEmbryoSelectedTissues_long[,"sampleType"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_PN","1-cell_embryo",humanEmbryoSelectedTissues_long[,"sampleType"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_ICM","blastocyst_ICM-polar_trophectoderm",
  humanEmbryoSelectedTissues_long[,"sampleType"])
humanEmbryoSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_TROPH","blastocyst_mural_trophectoderm",
  humanEmbryoSelectedTissues_long[,"sampleType"])

# order the sperm types (see fig 1 in Jan et al paper in /fh/fast/malik_h/grp/public_databases/NCBI/SRA/data/mammalian_expression_profiles/human/human_spermatogenesis)
humanEmbryoSelectedTissues_long[,"sampleType"] <- factor(
  humanEmbryoSelectedTissues_long[,"sampleType"],
  levels=c("1-cell_embryo", "2-4-8-cell_embryo", "morula", 
           "blastocyst_ICM-polar_trophectoderm", "blastocyst_mural_trophectoderm")) 



humanEmbryoSelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanEmbryoSelectedTissues_long[,"Gene"])
humanEmbryoSelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanEmbryoSelectedTissues_long[,"geneShort"])

humanEmbryoSelectedTissues_long[,"geneShort"] <- factor(
  humanEmbryoSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanEmbryoSelectedTissues_long[,"geneShort"]))   )

humanEmbryoSelectedTissues_long_summary <- humanEmbryoSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
humanEmbryo_plot <- ggplot(humanEmbryoSelectedTissues_long_summary, 
                                          aes(x=sampleType, y=my_median, fill=geneShort, 
                                              ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human embryo (Hendrickson)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))


########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/human_embryoPlots_2021_Mar24.pdf",sep=""), 
    width=7, height=7)
humanEmbryo_plot
dev.off()


######## human embryoOocyte
## (full version of the dataset used for human_embryo)

samplesOfInterest[["humanEmbryoOocyte"]][["human_embryoOocyte"]] <- 
  c("embryo_CL_rep1a", "embryo_CL_rep1b", "embryo_CL_rep2a", "embryo_CL_rep2b", 
    "embryo_ICM_rep1a", "embryo_ICM_rep1b", "embryo_ICM_rep2a", "embryo_ICM_rep2b", 
    "embryo_MOR_rep1a", "embryo_MOR_rep1b", "embryo_MOR_rep2a", "embryo_MOR_rep2b", 
    "embryo_PN_rep1a", "embryo_PN_rep1b", "embryo_PN_rep2a", "embryo_PN_rep2b",
    "embryo_TROPH_rep1a", "embryo_TROPH_rep1b", "embryo_TROPH_rep2a", "embryo_TROPH_rep2b",
    "immatureOocyte_GV_rep1a", "immatureOocyte_GV_rep1b", 
    "immatureOocyte_GV_rep2a", "immatureOocyte_GV_rep2b", 
    "immatureOocyte_MI_rep1a", "immatureOocyte_MI_rep1b", 
    "immatureOocyte_MI_rep2a", "immatureOocyte_MI_rep2b", 
    "matureOocyte_MII_rep1a", "matureOocyte_MII_rep1b", 
    "matureOocyte_MII_rep2a", "matureOocyte_MII_rep2b")

# get RPKMs for just the tissues of interest:
humanEmbryoOocyteSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                    species="humanEmbryoOocyte")

# select only some genes
humanEmbryoOocyteSelectedTissues <- humanEmbryoOocyteSelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanEmbryoOocyteSelectedTissues_long <- convertToLongFormat(humanEmbryoOocyteSelectedTissues)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoOocyteSelectedTissues_long[,"sample"])


# convert to long format
humanEmbryoOocyteSelectedTissues_long <- convertToLongFormat(humanEmbryoOocyteSelectedTissues)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoOocyteSelectedTissues_long[,"sample"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d[ab]$","",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_MOR","morula", humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_CL","2-4-8-cell_embryo",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_PN","1-cell_embryo",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_ICM","blastocyst_ICM-polar_trophectoderm",
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_TROPH","blastocyst_mural_trophectoderm",
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"])

unique(humanEmbryoOocyteSelectedTissues_long[,"sampleType"])

# order the embryo/oocyte types (see fig 1 in Hendrickson et al paper)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- factor(
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"],
  levels=c("immatureOocyte_GV", "immatureOocyte_MI", "matureOocyte_MII",
           "1-cell_embryo", "2-4-8-cell_embryo", "morula", 
           "blastocyst_ICM-polar_trophectoderm", "blastocyst_mural_trophectoderm")) 


humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanEmbryoOocyteSelectedTissues_long[,"Gene"])
humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanEmbryoOocyteSelectedTissues_long[,"geneShort"])

humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("H2Bopp8C","H2B.N",humanEmbryoOocyteSelectedTissues_long[,"geneShort"])


humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- factor(
  humanEmbryoOocyteSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanEmbryoOocyteSelectedTissues_long[,"geneShort"]))   )

humanEmbryoOocyteSelectedTissues_long_summary <- humanEmbryoOocyteSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
humanEmbryoOocyte_plot <- ggplot(humanEmbryoOocyteSelectedTissues_long_summary, 
                           aes(x=sampleType, y=my_median, fill=geneShort, 
                               ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human oocyte/embryo (Hendrickson)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))


########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/human_embryoOocytePlots_2021_Apr19.pdf",sep=""), 
    width=11, height=7)
ggarrange(humanEmbryoOocyte_plot, humanEmbryoOocyte_plot+coord_cartesian(ylim=c(0,15)),
          nrow=1)
dev.off()



######## human embryoOocyte
## (full version of the dataset used for human_embryo)

samplesOfInterest[["humanEmbryoOocyte"]][["human_embryoOocyte"]] <- 
  c("embryo_CL_rep1a", "embryo_CL_rep1b", "embryo_CL_rep2a", "embryo_CL_rep2b", 
    "embryo_ICM_rep1a", "embryo_ICM_rep1b", "embryo_ICM_rep2a", "embryo_ICM_rep2b", 
    "embryo_MOR_rep1a", "embryo_MOR_rep1b", "embryo_MOR_rep2a", "embryo_MOR_rep2b", 
    "embryo_PN_rep1a", "embryo_PN_rep1b", "embryo_PN_rep2a", "embryo_PN_rep2b",
    "embryo_TROPH_rep1a", "embryo_TROPH_rep1b", "embryo_TROPH_rep2a", "embryo_TROPH_rep2b",
    "immatureOocyte_GV_rep1a", "immatureOocyte_GV_rep1b", 
    "immatureOocyte_GV_rep2a", "immatureOocyte_GV_rep2b", 
    "immatureOocyte_MI_rep1a", "immatureOocyte_MI_rep1b", 
    "immatureOocyte_MI_rep2a", "immatureOocyte_MI_rep2b", 
    "matureOocyte_MII_rep1a", "matureOocyte_MII_rep1b", 
    "matureOocyte_MII_rep2a", "matureOocyte_MII_rep2b")

# get RPKMs for just the tissues of interest:
humanEmbryoOocyteSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                          species="humanEmbryoOocyte")

# select only some genes
humanEmbryoOocyteSelectedTissues <- humanEmbryoOocyteSelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanEmbryoOocyteSelectedTissues_long <- convertToLongFormat(humanEmbryoOocyteSelectedTissues)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoOocyteSelectedTissues_long[,"sample"])


# convert to long format
humanEmbryoOocyteSelectedTissues_long <- convertToLongFormat(humanEmbryoOocyteSelectedTissues)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanEmbryoOocyteSelectedTissues_long[,"sample"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d[ab]$","",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_MOR","morula", humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_CL","2-4-8-cell_embryo",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_PN","1-cell_embryo",humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_ICM","blastocyst_ICM-polar_trophectoderm",
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"])
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- gsub(
  "embryo_TROPH","blastocyst_mural_trophectoderm",
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"])

unique(humanEmbryoOocyteSelectedTissues_long[,"sampleType"])

# order the embryo/oocyte types (see fig 1 in Hendrickson et al paper)
humanEmbryoOocyteSelectedTissues_long[,"sampleType"] <- factor(
  humanEmbryoOocyteSelectedTissues_long[,"sampleType"],
  levels=c("immatureOocyte_GV", "immatureOocyte_MI", "matureOocyte_MII",
           "1-cell_embryo", "2-4-8-cell_embryo", "morula", 
           "blastocyst_ICM-polar_trophectoderm", "blastocyst_mural_trophectoderm")) 


humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanEmbryoOocyteSelectedTissues_long[,"Gene"])
humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanEmbryoOocyteSelectedTissues_long[,"geneShort"])

humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- gsub("H2Bopp8C","H2B.N",humanEmbryoOocyteSelectedTissues_long[,"geneShort"])


humanEmbryoOocyteSelectedTissues_long[,"geneShort"] <- factor(
  humanEmbryoOocyteSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanEmbryoOocyteSelectedTissues_long[,"geneShort"]))   )

humanEmbryoOocyteSelectedTissues_long_summary <- humanEmbryoOocyteSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
humanEmbryoOocyte_plot <- ggplot(humanEmbryoOocyteSelectedTissues_long_summary, 
                                 aes(x=sampleType, y=my_median, fill=geneShort, 
                                     ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human oocyte/embryo (Hendrickson)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))


########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/human_embryoOocytePlots_2021_Apr19.pdf",sep=""), 
    width=11, height=7)
ggarrange(humanEmbryoOocyte_plot, humanEmbryoOocyte_plot+coord_cartesian(ylim=c(0,15)),
          nrow=1)
dev.off()


######## human human_Xue_earlyEmbryo
## (full version of the dataset used for human_oocyte)


# gsub("_RPKM","",grep("RPKM", colnames(rpkms[["Hum_genes.blat_hg38.besthits.psl.bed.counts.human_Xue_earlyEmbryo.split.multicov.RPKMs.csv"]]) , value=TRUE))

# grep("RPKM", colnames(rpkms[["Hum_genes.blat_hg38.besthits.psl.bed.counts.human_Xue_earlyEmbryo.split.multicov.RPKMs.csv"]]) , value=TRUE)


samplesOfInterest[["humanXueOocyteEarlyEmbryo"]][["human_Xue_earlyEmbryo"]] <- 
  c("oocyte_rep1",             "oocyte_rep2",             "oocyte_rep3",            
    "pronucleus_rep1",         "pronucleus_rep2",         "pronucleus_rep3",     
    "zygote_rep1", "zygote_rep2",  
    "X2cell_blastomere_rep1",  "X2cell_blastomere_rep2",  "X2cell_blastomere_rep3",  
    "X4cell_blastomere_rep1",  "X4cell_blastomere_rep2",  "X4cell_blastomere_rep3", "X4cell_blastomere_rep4",  
    "X8cell_blastomere_rep1",  "X8cell_blastomere_rep2",  "X8cell_blastomere_rep3", "X8cell_blastomere_rep4",  
    "X8cell_blastomere_rep5",  "X8cell_blastomere_rep6",  "X8cell_blastomere_rep7", "X8cell_blastomere_rep8",  
    "X8cell_blastomere_rep9",  "X8cell_blastomere_rep10", "X8cell_blastomere_rep11",  
    "morula_rep1",             "morula_rep2",            "morula_rep3"         )

# get RPKMs for just the tissues of interest:
humanXueOocyteEarlyEmbryoSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                          species="humanXueOocyteEarlyEmbryo")
# select only some genes

humanXueOocyteEarlyEmbryoSelectedTissues <- humanXueOocyteEarlyEmbryoSelectedTissues[genesOfInterest[["human"]],]

# convert to long format
humanXueOocyteEarlyEmbryoSelectedTissues_long <- convertToLongFormat(humanXueOocyteEarlyEmbryoSelectedTissues)
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sample"])
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d+$","",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("^X","",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("pronucleus","pronuclear",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])


unique(humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])

# order the embryo/oocyte types (see fig 1 in Hendrickson et al paper)
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- factor(
  humanXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"],
  levels=c("oocyte", "pronuclear", "zygote", "2cell_blastomere", "4cell_blastomere", "8cell_blastomere", "morula")) 


humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("_Hum","",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"Gene"])
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("-pseudo","",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"])
humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("H2Bopp8C","H2B.N",humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"])

humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- factor(
  humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(humanXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"]))   )


humanXueOocyteEarlyEmbryoSelectedTissues_long_summary <- humanXueOocyteEarlyEmbryoSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
humanXueOocyteEarlyEmbryo_plot <- ggplot(humanXueOocyteEarlyEmbryoSelectedTissues_long_summary, 
                                 aes(x=sampleType, y=my_median, fill=geneShort, 
                                     ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="human oocyte/embryo (Xue)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
humanXueOocyteEarlyEmbryo_plot

########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/human_Xue_embryoOocytePlots_2021_Jun7.pdf",sep=""), 
    width=11, height=7)
ggarrange(humanXueOocyteEarlyEmbryo_plot, humanXueOocyteEarlyEmbryo_plot+coord_cartesian(ylim=c(0,15)),
          nrow=1)
dev.off()





######## mouse mouse_embryoOocyteXue
## (full version of the dataset used for mouse_embryo2)


# gsub("_RPKM","",grep("RPKM", colnames(rpkms[["Mouse_genes.blat_mouse_Dec2011.besthits.psl.bed.counts.mouse_embryoOocyteXue.split.multicov.RPKMs.csv"]]) , value=TRUE))

# Mouse_genes.blat_mouse_Dec2011.besthits.psl.bed.counts.mouse_embryoOocyteXue.split.multicov.RPKMs

samplesOfInterest[["mouseXueOocyteEarlyEmbryo"]][["mouse_embryoOocyteXue"]] <- 
  c("oocyte_rep1", "oocyte_rep2", 
    "pronucleus_rep1", "pronucleus_rep2", "pronucleus_rep3",   
    "X2cell_blastomere_rep1", "X2cell_blastomere_rep2", "X2cell_blastomere_rep3", 
    "X4cell_blastomere_rep1", "X4cell_blastomere_rep2", "X4cell_blastomere_rep3", 
    "X8cell_blastomere_rep1", "X8cell_blastomere_rep2", "X8cell_blastomere_rep3", 
    "morula_rep1", "morula_rep2", "morula_rep3")

# get RPKMs for just the tissues of interest:
mouseXueOocyteEarlyEmbryoSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                                  species="mouseXueOocyteEarlyEmbryo")
# select only some genes
mouseXueOocyteEarlyEmbryoSelectedTissues <- mouseXueOocyteEarlyEmbryoSelectedTissues[genesOfInterest[["mouse"]],]

# convert to long format
mouseXueOocyteEarlyEmbryoSelectedTissues_long <- convertToLongFormat(mouseXueOocyteEarlyEmbryoSelectedTissues)
mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sample"])
mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d+$","",mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])
mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("^X","",mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])
mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("pronucleus","pronuclear",mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])


unique(mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])

# order the embryo/oocyte types (see fig 1 in Hendrickson et al paper)
mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- factor(
  mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"],
  levels=c("oocyte", "pronuclear", "2cell_blastomere", "4cell_blastomere", "8cell_blastomere", "morula")) 


mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("_Mouse","",mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"Gene"])

mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- factor(
  mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(mouseXueOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"]))   )


mouseXueOocyteEarlyEmbryoSelectedTissues_long_summary <- mouseXueOocyteEarlyEmbryoSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
mouseXueOocyteEarlyEmbryo_plot <- ggplot(mouseXueOocyteEarlyEmbryoSelectedTissues_long_summary, 
                                         aes(x=sampleType, y=my_median, fill=geneShort, 
                                             ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="mouse oocyte/embryo (Xue)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
mouseXueOocyteEarlyEmbryo_plot

########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/mouse_Xue_embryoOocytePlots_2021_Jun7.pdf",sep=""), 
    width=7, height=7)
mouseXueOocyteEarlyEmbryo_plot
dev.off()





######## pig pig_bodyMap_batch1

# gsub("_RPKM","",grep("RPKM", colnames(rpkms[["Pig_genes.blat_susScr11.besthits.psl.bed.counts.pig_bodyMap_batch1.split.multicov.RPKMs.csv"]]) , value=TRUE))

samplesOfInterest[["pig_bodyMap_batch1"]][["pig_bodyMap_batch1"]] <- 
  c("adipose_1","adipose_2",
    "cerebellum_1","cerebellum_2","cerebellum_3",
    "cerebrum_1","cerebrum_2","cerebrum_3",
    "colon_1","colon_2","colon_3",
    "kidney_1","kidney_2","kidney_3",
    "liver_1","liver_2","liver_3",
    "lung_1","lung_2","lung_3",
    "muscle_biceps_1","muscle_biceps_2","muscle_biceps_3",
    "ovary_1","ovary_2","ovary_3",
    "smallIntestine_1","smallIntestine_2","smallIntestine_3",
    "spleen_1","spleen_2","spleen_3",
    "testis_1","testis_2","testis_3",
    "uterus_1","uterus_2","uterus_3")

# get RPKMs for just the tissues of interest:
pigBodymapBatch1SelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                                  species="pig_bodyMap_batch1")
# select only some genes
pigBodymapBatch1SelectedTissues <- pigBodymapBatch1SelectedTissues[genesOfInterest[["pig"]],]

# convert to long format
pigBodymapBatch1SelectedTissues_long <- convertToLongFormat(pigBodymapBatch1SelectedTissues)
pigBodymapBatch1SelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",pigBodymapBatch1SelectedTissues_long[,"sample"])
pigBodymapBatch1SelectedTissues_long[,"sampleType"] <- gsub("_\\d+$","",pigBodymapBatch1SelectedTissues_long[,"sampleType"])

unique(pigBodymapBatch1SelectedTissues_long[,"sampleType"])

# order the tissue types 
pigBodymapBatch1SelectedTissues_long[,"sampleType"] <- factor(
  pigBodymapBatch1SelectedTissues_long[,"sampleType"],
  levels=c("cerebellum","cerebrum",
           "adipose","muscle_biceps",
           "colon","kidney","liver","lung",
           "smallIntestine","spleen",
           "uterus","ovary","testis"))


pigBodymapBatch1SelectedTissues_long[,"geneShort"] <- gsub("_Pig","",pigBodymapBatch1SelectedTissues_long[,"Gene"])
pigBodymapBatch1SelectedTissues_long[,"geneShort"] <- gsub("-a$","",pigBodymapBatch1SelectedTissues_long[,"geneShort"])

pigBodymapBatch1SelectedTissues_long[,"geneShort"] <- factor(
  pigBodymapBatch1SelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(pigBodymapBatch1SelectedTissues_long[,"geneShort"]))   )


pigBodymapBatch1SelectedTissues_long_summary <- pigBodymapBatch1SelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
pigBodymapBatch1_plot <- ggplot(pigBodymapBatch1SelectedTissues_long_summary, 
                                         aes(x=sampleType, y=my_median, fill=geneShort, 
                                             ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="pig bodymap batch1") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
pigBodymapBatch1_plot

########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/pig_bodymapBatch1plots_2021_Jun7.pdf",sep=""), 
    width=11, height=7)
ggarrange(pigBodymapBatch1_plot, 
          pigBodymapBatch1_plot+coord_cartesian(ylim=c(0,10))+labs(title="pig bodymap batch1, zoom y-axis"),
          nrow=1)
dev.off()





######## pig pig_oocyteEarlyEmbryo

# gsub("_RPKM","",grep("RPKM", colnames(rpkms[["Pig_genes.blat_susScr11.besthits.psl.bed.counts.pig_oocyteEarlyEmbryo.split.multicov.RPKMs.csv"]]) , value=TRUE))

samplesOfInterest[["pig_oocyteEarlyEmbryo"]][["pig_oocyteEarlyEmbryo"]] <- 
  gsub("_RPKM","",grep("RPKM", colnames(rpkms[["Pig_genes.blat_susScr11.besthits.psl.bed.counts.pig_oocyteEarlyEmbryo.split.multicov.RPKMs.csv"]]) , value=TRUE))

# get RPKMs for just the tissues of interest:
pigOocyteEarlyEmbryoSelectedTissues <- collectRPKMsForGGplot(sampleList=samplesOfInterest, 
                                                         species="pig_oocyteEarlyEmbryo")
# select only some genes
pigOocyteEarlyEmbryoSelectedTissues <- pigOocyteEarlyEmbryoSelectedTissues[genesOfInterest[["pig"]],]

# convert to long format
pigOocyteEarlyEmbryoSelectedTissues_long <- convertToLongFormat(pigOocyteEarlyEmbryoSelectedTissues)
pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_RPKM","",pigOocyteEarlyEmbryoSelectedTissues_long[,"sample"])
pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("_rep\\d+$","",pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])


pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("^X","",pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])

pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- gsub("\\.","-",pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])


unique(pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"])

# order the tissue types 
pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"] <- factor(
  pigOocyteEarlyEmbryoSelectedTissues_long[,"sampleType"],
  levels=c("MII_oocyte",                
           "1-cell",                     
           "2-cell",                     
           "4-cell",                    
           "8-cell",                    
           "blastocyst-inner_cell_mass", 
           "blastocyst",                
           "blastocyst-trophoblast",   
           "morula" ))


pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("_Pig","",pigOocyteEarlyEmbryoSelectedTissues_long[,"Gene"])
pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- gsub("-a$","",pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"])

pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"] <- factor(
  pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"], 
  levels= intersect(preferredGeneOrder, 
                    unique(pigOocyteEarlyEmbryoSelectedTissues_long[,"geneShort"]))   )


pigOocyteEarlyEmbryoSelectedTissues_long_summary <- pigOocyteEarlyEmbryoSelectedTissues_long %>% 
  group_by(sampleType,geneShort) %>% 
  summarise(my_median=median(RPKM), my_mad=mad(RPKM))

### using summarized data:
pigOocyteEarlyEmbryo_plot <- ggplot(pigOocyteEarlyEmbryoSelectedTissues_long_summary, 
                                aes(x=sampleType, y=my_median, fill=geneShort, 
                                    ymin=my_median-my_mad, ymax=my_median+my_mad)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(stat="identity", position=position_dodge(0.9), width=0.5, color="gray") + 
  genesFillScale +
  guides(fill=guide_legend(title="gene")) +
  labs(x="", y="RPKM", title="pig oocyte/early embryo (Kong et al 2020)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5))
pigOocyteEarlyEmbryo_plot

########## save plot and zoomed version:
pdf(file=paste(plotsDirPretty,"/pig_oocyteEarlyEmbryoPlots_2021_Jun16.pdf",sep=""), 
    width=7, height=7)
pigOocyteEarlyEmbryo_plot
dev.off()


