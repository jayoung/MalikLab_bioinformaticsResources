
### set up a list object that has the types of sequence we see at each "site"
bitsToCombine <- list()
bitsToCombine[[1]] <- "leftFlankingSeq"
bitsToCombine[[2]] <- c("A","B")
bitsToCombine[[3]] <- c("C","D")
bitsToCombine[[4]] <- c("E","F")
bitsToCombine[[5]] <- "middlePadding"
bitsToCombine[[6]] <- c("G","H")
bitsToCombine[[7]] <- c("I","J")
bitsToCombine[[8]] <- "rightFlankingSeq"


## go through the list

# initialize an empty object
sequenceCombos <- ""
expectedNumCombos <- 1

for(i in 1:length(bitsToCombine)) {
    thisBitOptions <- bitsToCombine[[i]]
    howManyOptions <- length(thisBitOptions)
    expectedNumCombos <- expectedNumCombos * howManyOptions
    howManyCombosExistAlready <- length(sequenceCombos)
    cat("round",i,"repeating the initial list",howManyOptions,"times and adding these choices",thisBitOptions,"\n")
    
    # multiply the existing combos
    sequenceCombos <- rep(sequenceCombos, howManyOptions)
    
    # repeat the current site options enough times to add it to each one of the strings we are building up
    bitsToAdd <- rep(thisBitOptions, each=howManyCombosExistAlready)

    # check we didn't mess up
    if (length(sequenceCombos) != length(bitsToAdd)) {
        stop("ERROR - something is wrong with the code\n")
    }
    sequenceCombos <- paste(sequenceCombos, bitsToAdd, sep="_")
}

# take a look
sequenceCombos

length(unique(sequenceCombos))

# expected length
expectedNumCombos
