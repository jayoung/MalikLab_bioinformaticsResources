#!/usr/local/bin/python3
import sys
import Bio

# import all the NCBI web methods
from Bio.Blast import NCBIWWW

fasta_string = open("humanCENPAcds.fa").read()

## do the blast (result_handle type=_io.StringIO)
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

# (default output format is XML)
help(NCBIWWW.qblast)

## get the result (blast_record type=Bio.Blast.Record.Blast)
from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

len(blast_record.alignments) 
type(blast_record.alignments) # list

# parse through the hits (first 3 only, all of which satisfy the e-value threshold)
E_VALUE_THRESHOLD=0.01

for alignment in blast_record.alignments[0:3]:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESHOLD:
            print ("\n### alignment")
            print("hit name: ", alignment.title)
            print("length:   ", alignment.length)
            print("e-value:  ", hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

