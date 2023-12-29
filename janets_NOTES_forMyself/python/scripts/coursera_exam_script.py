#!/usr/local/bin/python3
import sys
import getopt
import Bio
import re

def usage():
    print ("""
coursera_exam_script.py : reads fasta file, and answers the following questions:
  1. how many sequences are in the file
           
  2a. what are the lengths of each sequence (we will output a tab-delim text file of all lengths)
  2b. What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers? 
           
  3a. identify all ORFs present in each sequence of the FASTA file (only the forward strand, i.e. frames 1, 2, or 3) 
  3b. what is the length of the longest ORF in the file? 
  3c. What is the identifier of the sequence containing the longest ORF? 
  3d. For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? 
  3e. What is the starting position of the longest ORF in the sequence that contains it? 
           
  4a. Given a length n, identify all repeats of length n in all sequences in the FASTA file. 
  4b. Determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.

coursera_exam_script.py [-h] [-k <kmerSize>] <filename>
  -h             print this message
  -s <seqname>   name of the sequence we want to get the longest ORF for
  -k <kmerSize>  kmer size for repeat analysis <kmerSize>
  <filename>     filename in fasta format
    """)

### get options
# o = optional args (a list)
# a = required args (a list)
# the l:h says there could be two possible options, h and l. The colon after the l means we expect a value after the l
o, a = getopt.getopt(sys.argv[1:], 's:k:h')

opts = {}
# set defaults:
kmer=0
seqname=""

# create dictionary for the optional args
for k,v in o:
    opts[k]=v

if '-h' in opts.keys():
    usage(); sys.exit()

if '-k' in opts.keys():
    opts['-k']=int(opts['-k'])
    kmer=opts['-k']

if '-s' in opts.keys():
    seqname=opts['-s']

if len(a) < 1:
    usage(); sys.exit("\nERROR - fasta file name is missing\n\n")
mySeqFile=a[0]

## make sure the file exists
try:
    f = open(mySeqFile, "r")
except IOError: 
    print("\n\nERROR: file %s does not exist\n\n" % mySeqFile)

f.close()

## figure out filename for the output file we'll use for seq lengths, and open it for writing
lenFileName = re.sub("\.fa$", "", mySeqFile)
lenFileName = re.sub("\.fasta$", "", lenFileName)
lenFileName = lenFileName + ".lengths.tsv"
lenFile = open(lenFileName, "w")
lenFile.write("SeqID\tLength\n") 


## figure out filename for the output file we'll use for longest ORFs and open it for writing
orfFileName = re.sub("\.fa$", "", mySeqFile)
orfFileName = re.sub("\.fasta$", "", orfFileName)
orfFileName = orfFileName + ".ORFs.tsv"
orfFile = open(orfFileName, "w")
orfFile.write("SeqID\tLength\tORFlen\tFrame\tStart\tStop\n") 

## define a function to check for longest ORF
def getAllORFs(inputSeq, frame):
    frameStartPos = frame-1
    myseq=inputSeq[frameStartPos:]
    startPositions=[]
    stopPositions=[]
    allORFs=[]
    withinAnOrf=0
    all_start_codons=["atg","ATG"]
    all_stop_codons=["taa","tga","tag","TAA","TGA","TAG"]
    for i in range(0,len(myseq),3):
        thisCodon=str(myseq[i:(i+3)])
        ## check for start codon.  If we're already within an ORF we don't do anything
        if thisCodon in all_start_codons:
            if withinAnOrf==0:
                startPositions.append(i)
                withinAnOrf=1
        ## check for stop codon.  If we're NOT within an ORF we don't do anything
        if thisCodon in all_stop_codons:
            if withinAnOrf==1:
                stopPositions.append(i)
                withinAnOrf=0
    ## finished looping through codons
    ## did we find at least one ORF?
    # if so, we process each ORF and only retain the longest
    if len(startPositions)>0 and len(stopPositions)>0:
        for i in range(len(startPositions)):
            thisStart = startPositions[i]
            # this if statement means we ignore starts that don't have a matching stop
            if i<len(stopPositions):
                thisStart = startPositions[i]
                thisStop = stopPositions[i]
                ORFlen = thisStop+3-thisStart
                ORFseq = myseq[thisStart:(thisStop+3)]
                ORFobject = {"len":ORFlen, "seq":ORFseq, "start":(thisStart+1+frameStartPos), "stop":(thisStop+3+frameStartPos), "frame":frame}
                allORFs.append(ORFobject)
    return(allORFs)

## read in seqs
from Bio import SeqIO
seqdict = {}
longestSeqLen = 0
longestSeqNames = []
allKmers={} # key will be the kmer seq, value will be the count
for seq_record in SeqIO.parse(mySeqFile, "fasta"):
    seqid = seq_record.id
    # store sequence in the dictionary
    seqdict[seqid]=seq_record
    # write seq length to file
    lenString = seqid + "\t" + str(len(seq_record)) + "\n"
    lenFile.write(lenString) 
    # check whether it's longer than longest thing we've seen yet
    if len(seq_record) > longestSeqLen:
        longestSeqLen = len(seq_record)
        longestSeqNames = [seqid]
    # check whether it's the same length as the longest thing we've seen yet
    elif len(seq_record) == longestSeqLen:
        longestSeqNames.append(seqid)
    ##### kmers - count them!
    if kmer > 1:
        lastKmerStartPos = len(seq_record) +1 - kmer
        for i in range(lastKmerStartPos):
            thisKmer=str(seq_record.seq[i:(i+kmer)])
            if thisKmer in allKmers.keys():
                allKmers[thisKmer]=allKmers[thisKmer]+1
            else:
                allKmers[thisKmer]=1
    ##### ORFs
    # check for ORFs in each of the forward strand reading frames
    allORFs=[]
    allORFs.append( getAllORFs(seq_record.seq, frame=1) )
    allORFs.append( getAllORFs(seq_record.seq, frame=2) )
    allORFs.append( getAllORFs(seq_record.seq, frame=3) )
    # go through all ORFs to get the longest one(s)
    longestORFlen=0
    longestORFobjects=[]
    for i in range(len(allORFs)):
        ORFsThisFrame=allORFs[i]
        for j in range(len(ORFsThisFrame)):
            thisORF=ORFsThisFrame[j]
            # if it is LONGER:
            if thisORF["len"] > longestORFlen:
                longestORFlen=thisORF["len"]
                longestORFobjects=[thisORF]
            elif thisORF["len"] == longestORFlen:
                longestORFobjects.append(thisORF)
    if longestORFlen>0:
        # it's unlikely but it is possible there's a tie for longest ORF
        numLongestORFs=len(longestORFobjects)
        # print("\nin seq %s (length %d) there was %d longest ORF:" % (seqid, len(seq_record), len(longestORFobjects) ))
        for i in range(len(longestORFobjects)):
            # print("frame %d, ORF len %d, pos %d - %d" % (longestORFobjects[i]["frame"], longestORFobjects[i]["len"], longestORFobjects[i]["start"],longestORFobjects[i]["stop"] ))
            orfOutputString = "\t".join([seqid, str(len(seq_record)), str(longestORFobjects[i]["len"]), str(longestORFobjects[i]["frame"]), str(longestORFobjects[i]["start"]), str(longestORFobjects[i]["stop"])])
            orfOutputString=orfOutputString+"\n"
            orfFile.write(orfOutputString) 
            # print("ORF seq %s" % longestORFobjects[i]["seq"])
            # translation = longestORFobjects[i]["seq"].translate()
            # print("pep seq %s" % translation)
    else:
        print("in seq %s  (length %d) frame1 had no ORFs" % (seqid, len(seq_record)))

## How many seqs are in the file
numseqs=len(seqdict)
print ("\nThere are %d seqs in file %s\n" % (numseqs,mySeqFile))

## show longest seq(s)
longestSeqNamesString = ", ".join(longestSeqNames)
print("The longest seq(s) were %d bp long and is/are called %s\n" % (longestSeqLen,longestSeqNamesString))


## show repeats
foundAny=0
mostFrequentRepeatCount=0
mostFrequentRepeats=[]
if kmer > 1:
    # go through each kmer
    for kmerSeq in allKmers.keys():
        kmerCount=allKmers[kmerSeq]
        if kmerCount>1:
            print ("found kmer %s %d times" % (kmerSeq,kmerCount))
            foundAny = foundAny+1
            if kmerCount>mostFrequentRepeatCount:
                mostFrequentRepeatCount=kmerCount
                mostFrequentRepeats=[kmerSeq]
            elif kmerCount==mostFrequentRepeatCount:
                mostFrequentRepeats.append(kmerSeq)
    if foundAny==0:
        print("did not find any repeated %d-mers" % kmer)
    else:
        freqRepeatString = "\n".join(mostFrequentRepeats)
        print("\nThe most frequent repeat(s) were found %d times and they are:\n%s" % (mostFrequentRepeatCount,freqRepeatString))


