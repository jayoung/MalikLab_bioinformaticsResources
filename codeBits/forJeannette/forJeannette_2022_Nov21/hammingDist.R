library(Biostrings)

seqs1 <- c("CNIVPLN",
          "CNIAPLN",
          "CNIVPWN",
          "CNIVPWN")
seqs2 <- AAStringSet(seqs1)

##### get all pairwise hamming distances
hammingMatrix <- sapply( seqs2, function(x) {
    neditAt(x, seqs2, at=1)
} )
rownames(hammingMatrix) <- seqs1
colnames(hammingMatrix) <- seqs1
