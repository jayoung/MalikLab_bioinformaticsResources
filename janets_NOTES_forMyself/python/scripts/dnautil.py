#!/usr/local/bin/python3
"""
a few functions to do stuff to DNA seqs
"""

def gc(dna) :
    """
    compute GC percentage
    maybe I add more info here
    """ 
    nbases=dna.count("n") + dna.count("N")
    gcpercent=float( dna.count("c") + dna.count("C") + dna.count("g") + dna.count("G") )*100.0 / (len(dna)-nbases)   
    return gcpercent

def check_for_stops(dna, frame=0) :
    "a function to check a DNA sequence for in-frame stop codons"
    any_stop_codon = ["taa","tga","tag"]
    has_stop = False
    for startPos in range(frame,len(dna),3) :
        thisCodon = dna[ startPos : (startPos+3) ].lower()
        if thisCodon in any_stop_codon : 
            has_stop = True
            break
    return(has_stop)

def revcomp (seq) :
    """Return reverse-complement of a DNA sequence"""
    # first we reverse the string:
    seq = seq[::-1]  # same as [len(seq):0:-1]
    # then we complement
    baseComplementDict = {"A":"T", "C":"G", "G":"C", "T":"A", "a":"t", "c":"g", "g":"c", "t":"a", "n":"n", "N":"N"}
    # coerce to list: this splits up the letters
    letters = list(seq)
    # the following is a "list comprehension"
    letters = [baseComplementDict[eachBase] for eachBase in letters]
    # concatenate them back together
    newSeq = "".join(letters)
    return newSeq 
