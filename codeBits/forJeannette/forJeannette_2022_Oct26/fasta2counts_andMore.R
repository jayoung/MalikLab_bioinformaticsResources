## Question: 
# 
# (2) After all my pre-processing (trimming + read merging + quality filtering), I used FastX-collapser to collapse all my reads down to a fasta file with the # of counts for unique sequence (in an effort to avoid loading huge files into R). This gave me an awkward fasta file that looks like this:
# >1_counts
# seq
# >2_counts
# seq 
# And I’m not quite sure how to extract that into a table with a seq and counts column. Perhaps there’s a better way to do this in Terminal (with sort uniq?); fastX-collapser was just easy and quick.


### try it!

# Biostrings package has the functions we need

library(Biostrings)

# I made an example fasta file:
fastaFile <- "forJeannette_2022_Oct26/exampleFastaWithCounts.fa"

# read in fasta file
fa <- readDNAStringSet(fastaFile)
fa

# we can access the sequence and the names
dna <- as.character(fa)
dna

fa_headers <- names(fa)
fa_headers

# we can split the ids and get the component parts
fa_headers_split <- strsplit(fa_headers, "_")

seq_ids <- sapply(fa_headers_split, "[[", 1)
seq_counts <- as.integer(sapply(fa_headers_split, "[[", 2))

# putting all that code together inside a function, and producing a data frame at the end:
fasta2counts <- function(file) {
    fa <- readDNAStringSet(file)
    fa_headers_split <- strsplit(names(fa), "_")
    output <- data.frame(id=sapply(fa_headers_split, "[[", 1),
                         counts=as.integer(sapply(fa_headers_split, "[[", 2)),
                         sequence=as.character(fa))
    return(output)
}
# use that function
countsTable <- fasta2counts(fastaFile)


# might want to save it as tab-delimited text
write.table(countsTable, 
            file="forJeannette_2022_Oct26/exampleFastaWithCounts.countsTable.tsv", 
            row.names=FALSE, quote=FALSE,
            sep="\t")



####### translation

AAs <- translate( DNAStringSet(countsTable$sequence), 
                  no.init.codon=TRUE )
AAs

# note:  I use the "no.init.codon=TRUE option" because Biostrings' translate function has a weird (to me) default setting.  It tries hard to translate the first codon in whatever sequences you give it as an initiation codon, so if the first three letters are ATG, CTG or TTG, you will get M when you translate.  (CTG and TTG are alternative start codons for a minority of genes). You can turn that behavior off using the “no.init.codon=TRUE option”

# need to convert back to character if we want to add it to our data.frame
countsTable$AA <- as.character(AAs)

######### define the distance from WT in terms of nt changes and amino acid changes

# first we define the WT seqs. For this example let's pretend it's the first one in our counts table. I'm sure it's not in the real example.
WT_dna <- DNAString(countsTable[1,"sequence"])
WT_aa <- translate(WT_dna, no.init.codon=TRUE)

# IF there are no indels (the DMS only library) we don't need to align seqs before calculating distance

# Biostrings' neditAt function calculates hamming distance
NTdists <- neditAt(WT_dna, 
                   DNAStringSet(countsTable$sequence))
NTdists

# NT dists is in an annoying format but that easy to deal with:
t(NTdists)

## add NT and AA dist to table:
countsTable$NTdistFromWT <- t(NTdists)
    
countsTable$AAdistFromWT <- t( neditAt(WT_aa, 
                                    AAStringSet(countsTable$AA)) )



####### get position of the first amino acid change

# there's probably a way to do this in Biostrings, but it's more obvious to me how to do it in plain R:

## we can split seqs into individual letters like this

# the strsplit function returns a LIST object (each item is a character vector)
AAs_split <- strsplit(as.character(AAs), "") 
AAs_split

# for WT, there's only 1 seq, so we can take the first list element using [[1]]
WT_aa_split <- strsplit(as.character(WT_aa), "")[[1]]
WT_aa_split

# now, for an example item (the second, as the first IS wt) in AAs_split, we'll see where the first mismatch is

# get all mismatches
test <- WT_aa_split != AAs_split[[2]]

# which() shows positions of TRUE items:
which(test)

# and this is the first mismatch
which(test)[1]

### put that together in a function. Operates on character objects. Should work for NT or AA seqs
getFirstChangePosSplitSeqLists <- function(wtSeq, querySeqs) {
    # split the seqs:
    wtSeq_split <- strsplit(as.character(wtSeq), "")[[1]]
    querySeqs_split <- strsplit(as.character(querySeqs), "")
    # for EACH item in querySeqs_split list, we get the first change pos
    firstChangePositions <- sapply( querySeqs_split, function(thisSplitSeq)  {
        test <- wtSeq_split != thisSplitSeq
        firstChange <- which(test)[1]
        return(firstChange)
    })
    return(firstChangePositions)
}
countsTable$firstNTchangePos <-getFirstChangePosSplitSeqLists( 
    as.character(WT_dna),
    countsTable$sequence )

countsTable$firstAAchangePos <-getFirstChangePosSplitSeqLists( 
    as.character(WT_aa),
    as.character(AAs) )




