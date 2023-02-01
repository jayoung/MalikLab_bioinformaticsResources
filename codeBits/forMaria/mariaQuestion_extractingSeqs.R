library(Biostrings)
library(GenomicRanges)
library(BSgenome) # needed for the getSeq function

## make example transcript sequences
transcripts_plain <- c("ACGTGATCGTGGACTGTAGCTGG",
                       "AGCTGACTGGTGTTGTTTTTT",
                       "TGTTGATGGTCGTGATGTTTAAAAA")
names(transcripts_plain) <- c("seq1", "seq2", "seq3")

## convert sequences to a DNAStringSet object (Biostrings package)
transcripts <- DNAStringSet(transcripts_plain)

## make example regions
desiredRegions <- GRanges(seqnames=c("seq1","seq2","seq3"),
                          ranges=IRanges(start=c(4,4,4), end=c(10,10,10)))

## use getSeq function
getSeq(transcripts, desiredRegions)
