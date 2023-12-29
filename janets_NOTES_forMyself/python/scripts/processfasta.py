#!/usr/local/bin/python3
import sys
import getopt

def usage():
    print ("""
processfasta.py : reads fasta file, builds dictionary

processfasta.py [-h] [-l <length>] <filename>
  -h           print this message
  -l <length>  remove seqs with length > than <length>
  <filename>   filename in fasta format
    """)

### use getopt!
# o = optional args (a list)
# a = required args (a list)
# the l:h says there could be two possible options, h and l. The colon after the l means we expect a value after the l
o, a = getopt.getopt(sys.argv[1:], 'l:h')

opts = {}
seqlen=0 # default

# create dictionary for the optional args
for k,v in o:
    opts[k]=v

if '-h' in opts.keys():
    usage(); sys.exit()

if len(a) < 1:
    usage(); sys.exit("\nERROR - fasta file name is missing\n\n")

if '-l' in opts.keys():
    opts['-l']=int(opts['-l'])
    if opts['-l']<0:
         sys.exit("\nERROR - can't have -l <0\n\n")
    seqlen=opts['-l']

mySeqFile=a[0]

try:
    f = open(mySeqFile, "r")
except IOError: 
    print("\n\nERROR: file %s does not exist\n\n" % mySeqFile)

# read each line - is it a header?
mydict = {}
for line in f: 
    line=line.rstrip() # "return strip"
    ## there's also a startswith method we could have used instead of [0]
    if line[0] == ">":
        # if so get name, make new dictionary entry
        id=line.split()[0][1:]  # the 1: is to get rid of the >
        mydict[id]=""
    else:
        # if not, append seq to dicitionary value
        mydict[id]=mydict[id]+line

# close file
f.close()

# or, fancier:
for name,seq in mydict.items():
    if len(seq)<=seqlen:
        print(name,seq)
    else:
        print("seq",name,"is too long")

