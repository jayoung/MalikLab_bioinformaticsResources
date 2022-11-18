######### plot preferences using DiffLogo (differences, not ratios) 
# see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4650857/

library(DiffLogo)
library(motifStack)

countTableFiles <- c("forJeannette_2022_Nov18/rossana_countTable_library.csv",
                     "forJeannette_2022_Nov18/rossana_countTable_super-restrictor.csv")
countTables <- lapply(countTableFiles, function(x) {
    y <- read.table(x, sep=",", row.names = 1, header=TRUE)
    y <- t(y) # transpose (sites should be columns not rows)
    return(y)
})
countTables

## convert to frequency tables, not raw counts
# these are matrices
freqTables <- lapply(countTables, function(x) {
    z <- apply(x, 2, function(y) { y/sum(y) })
    return(z)
} )


### try diffLogo.  The error message is TOTALLY misleading - I'll walk through why that is.
diffLogoFromPwm(freqTables[[1]], freqTables[[2]])
# [1] "pwm must be of class matrix or data.frame. Trying to convert"
# [1] "pwm must be of class matrix or data.frame. Trying to convert"
# Error in preconditionPWM(pwm1) : Columns of PWM must add up to 1.0

### convince myself that columns DO add up to 1:
lapply(freqTables, colSums)

### turns out the problem is not about the column sums, it's that the default is that it expects a DNA matrix, i.e. only 4 rows, so it's confused by us giving it a matrix with 20 rows.
# we use the 'alphabet' argument to tell it how many rows there should be. This was very non-obvious from the help.  
diffLogoFromPwm(freqTables[[1]], 
                freqTables[[2]], 
                alphabet=ASN)

### some more notes on alphabet - might need it to handle the +/- 2 aa indel

# ASN alphabet has the usual 20 amino acids and an associated color scheme. See what that looks like:
ASN

# FULL_ALPHABET is the 26 letters of the English language alphabet, also with some colors.  I did use it as the basis for messing around with the colors in the example code below...
FULL_ALPHABET

#### here's some more code I used because I wanted to mess with the colors.   You might want something like this to help you with the indel

### define some utility functions and setup to use DiffLogo:

### make a color scheme that mimics the "chemistry" color scheme from motifStack
# see https://github.com/jianhong/motifStack/blob/5aa80388b44bc8f93738315310cfb56c3495130c/R/publicUtilities.R 
# with this code to figure out how to convert #800080 colors into text colors
# plot(1:10,1:10, pch=19, cex=2, col="#800080")
# points(1:10,0:9, pch=19, cex=3, col="magenta4")

## changeColors: a little function to change color for a single component of the alphabet
changeColors <- function(myAlphabet, myAA, myCol) {
    whichChar <- which(myAlphabet$chars %in% myAA)
    tempColors <- myAlphabet$cols
    tempColors[whichChar] <- myCol
    myAlphabet$cols <- tempColors
    myAlphabet
}
## reportColors: a little function just to help me look at the color assignments I made
reportColors <- function(myAlphabet) {
    x <- data.frame(aa=myAlphabet$chars, color=myAlphabet$cols)
    print (x)
    y <- table(x[,"color"])
    print (y)
}

## now we take FULL_ALPHABET and mess with the colors
FULL_ALPHABET_JYchemistryColors <- FULL_ALPHABET
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("A","F","I","L","M","P","V","W"), "black")
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("C","G","S","T","Y"), "forest green")
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("D","E"), "red3")
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("H","K","R"), "blue3")
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("N","Q"), "magenta4")
# and pink for the amino acids that don't exist
FULL_ALPHABET_JYchemistryColors <- changeColors(FULL_ALPHABET_JYchemistryColors, c("B","J","O","U","X","Z"), "pink")
reportColors(FULL_ALPHABET_JYchemistryColors)

## reformatPFMincludeFullAlphabet: a function to reformat the PFMs so that they have rows for ALL the amino acids mentioned in the alphabet (the alphabet includes rows for amino acids that don't exist in my counts/freqs tables, and that causes trouble)
reformatPFMincludeFullAlphabet <- function(myPWM) {
    newDat <- data.frame(row.names=FULL_ALPHABET$chars)
    numCols <- dim(myPWM)[2]
    for(tempColIndex in 1:numCols ) {
        newDat[ rownames(myPWM),tempColIndex] <- myPWM[,tempColIndex]
        newDat[which(is.na(newDat[,tempColIndex])),tempColIndex] <- 0
    }
    newDat
}


freqTables_FULL_ALPHABET <- lapply( freqTables, reformatPFMincludeFullAlphabet)

### difference logo
diffLogoFromPwm( freqTables_FULL_ALPHABET[[1]], 
                 freqTables_FULL_ALPHABET[[2]], 
                 alphabet = FULL_ALPHABET_JYchemistryColors)  


