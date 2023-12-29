# Coursera: Python for genomic data science

"Python for genomic data science" [course](https://www.coursera.org/learn/python-genomics/home/week/1).

## Course contents
- module1:  
  - lecture1: Overview of Python   
  - lecture2: First Steps Toward Programming  
- module2:
  - lecture3: Data Structures
  - lecture4: Ifs and Loops 
- module3:
  - lecture5: Functions
  - lecture6: Modules and Packages
- module4:
  - lecture7: Communicating with the Outside
  - lecture8: Biopython
- exam:
  - see `coursera_exam_notes.md` and `coursera_exam_script.py` (and input file `dna.example.fasta`)


## Suggested resources
'Beginning python for bioinformatics' by Patrick O'Brien

learnpython.org

thinkpython online book - think like a computer scientist
https://docs.python.org/3.9/tutorial/

Biopython [Tutorial/cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html)

Biopython [FAQs](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec6)

## My notes

python2 versus python3 - There's some weirdness about python2:  12/5=2, because it confuses integer and real numbers. 12 is integer, so the output to the function is also an integer. in contrast 12.0/5=2.4 and float(12)/5=2.4 

Use python3!

Comments:
```
"""
comments across several lines can 
be shown using triple-quotes
"""
```
For shorter comments we can use the # char (anywhere in a line)

### Data types
- **numbers**
   - **integer** 
   - **float** ("float" = "real")  
   - complex/imaginary numbers (also available, use j) (not often used in biology

- **strings**
    - single quote or double-quote: double quotes better because you might want to use an apostrophe e.g. "it's"
    - use \ to escape quotes or whatever else
    - use triple quotes for strings that span multiple lines (it interprets newlines and converts them to \n
    - common escaped things: \n \t \\ \"
    
- **list**  
ordered set of values, can be mixed types

- **tuple**  
similar to lists but immutable

- **sets**  
an unordered collection with no duplicates (can use union, intersection, difference)

- **dictionaries**  
key-value pairs

Functions versus methods:

e.g. `len` is a function, and we call it like this:
```
dna1 = "atgtagctgcgatcgtgatcg"
len(dna1)

len
# <built-in function len>
```
but `upper` is a method, and we call it like this:
```
dna1.upper()

dna1.upper
# <built-in method upper of str object at 0x10a69ed50>
```

## Some code bits

Start up an interactive python session
```
python3
```

My first command:
```
print("hello")
```

More commands:
```
# numbers:
5+5

# "**" means "to the power of"
10**2  

# "//" means "floor division" (discard the fraction)
17//3

# "%" means modulo (return the remainder after division)
17 % 3 

## 'type' function is like class in R
type(5)
type(5.0)
type("hello")
```

### Using quotes in python scripts:
```
# simply putting this on the command line doesn't print it nicely:
"""\
>dna1
agtgatagtgtgcgcgtgcggtag
agggtgggtagtc
>dna1
agtgatagtgtgcgcgtgcggtag
agggtgggtagtc
"""

# but using the print function it does print better:
print("""\
>dna1
agtgatagtgtgcgcgtgcggtag
agggtgggtagtc
>dna1
agtgatagtgtgcgcgtgcggtag
agggtgggtagtc
""")
```

### "strong operaters"

+, *, in, not in
```
# + (concatenate)
"aaa" + "ggg"
# * (replicate)
"atg" * 3
# in / not in
"el" in "hello" # (returns True, of class 'bool')
"el" not in "hello" # False
"nel" in "hello" # True
```

### Variables 

In variable names, numbers are OK if they're NOT in the first position. No special chars except _
```
codon = "atg"
dna = "atgtagctggatcgtgatcg"
```

### Substrings 

0 is the first position

```
codon[0] # first char
codon[1] # second char
codon[1:2] # end position is always excluded (just like bed format) (second and third)
codon[-1] # count backwards 1 from end of string
codon[1:] # from second char to end
codon[:2] # from start to pos 3
codon[:3]
```

### Functions

`len` (length of a string)
```
help(len) # within python
pydoc len # from a shell (but not on my mac)
```


### Object-oriented programming

e.g. there are some functions that only objects of a specific type can do. e.g. `count` is a function for strings
```
dna.count("c")
    # the period means: "ask object dna to do the count function"
dna.count("gc")
dna.upper()
dna.lower()
```
find function: 
```
# find the first occurrence (gives 0-based position)
dna.find("g")
# second arg means start looking from a certain position onwards
dna.find("g",3)
# find starting at the end of the string
dna.rfind("g")
```

show all methods available for string objects: `help(str)`

More string functions:
```
# are all chars in dna in uppercase?
dna.isupper()
# replace all instances:
dna.replace("a","A")
```


## Writing scripts!

Simple script to get GC content - see `gc.py`


## `input` function 

(=`raw_input` in python2)

`input()` always captures the input as a string, so often we want to convert that to something else, e.g. 
```
my_number = input("Please type in a number:")
my_number = int(my_number)
```


## conversion functions
```
int()
int(x, [,base])
float(x)
complex(real[,imag])
str(x) # converts integer to a string
chr(x) # converts integer to a character (i.e. ascii)

str(65) # result = "65"
chr(65) # result = "A"
```

## print with formatting

```
print("my message") # paretheses not needed in python2 but they are in python3
```
This primntes a number nicely:
```
print("gcPercent is %5.1f %%" % gcPercent)
```
`%` says a formatting command follows  
`5` is for total number of digits  
`1` is for num decimal places  
`f` is for floating point  

```
print ("%d" % 10.6)   # %d means transform to an integer
print ("%5d" % 10.6)  # result = "   10" - the 5 results in some padding using spaces to get total of 5 digits

print ("%o" % 10) # convert to octal. result=12 (fails on a float, e.g. 10.6)
print ("%x" % 10) # convert to hexadecimal result=a (fails on a float)

print ("%e" % 10) # show in scientific notation (result = 1.060000e+01)

print ("%s" % dna) # show in scientific notation
```

## lists
```
gene_expression = ['gene', 5.16e-08, 0.000138511, 7.33e-08]
gene_expression[0] = "Lif"
len(gene_expression)
```

lists are mutable data types but strings are immutable (so the following reassignment `dna[0] = "n"` would NOT work)

`gene_expression[:]` is a special 'slice' operation returns the whole list  
`gene_expression[:] = []`   # empty list

concatenate lists (without modifying)
```
gene_expression + ['newElement1','newElement2']
```
del - remove elements (destructively) (shifts remaining elements to the left).  Doesn't look like we're changing gene_expression here but we are:
```
del gene_expression[1] 
```


list-related methods: `extend`, `count`:
```
gene_expression.extend(['newElement1','newElement2'])
gene_expression.count("Lif")  # exact matches to complete list elements (e.g. 'Lif' is there 1 time, 'Li' is there 0 times)

print(gene_expression.count("Lif"), gene_expression.count("gene"))
```

`reverse()` - modifies gene_expression by reversing it
```
gene_expression.reverse()
```

help(list)

append and pop:
```
append
pop

stack = ["a","b","c","d"]

stack.append("e")

# appending a list results in a NESTED list in the last element, so it's different from extend/concatenation
stack.append(["f","g"])

# pop removes the last list element
lastThing = stack.pop()
# remove but don't save the last list element
stack.pop()
```

sorting lists
```
mylist = [3,31,123,1,5]

# sorts and outputs to screen but doesn't change mylist
sorted(mylist)

# DOES change mylist:
mylist.sort()

# sort can work on non-numerical data
mylist = ['a','g','a','c','a']
sorted(mylist)

# sort fails on lists of mixed type
mylist = ['g','a','c','a',10,1,5,50]
sorted(mylist)
```

## tuples

Like lists, but immutable (?)
```
t=1,2,3
t=(1,2,3) # same as above
```
can use most of same functions on tuple as we could on list as long, as they don't change anything


## sets 
A set is an unordered collection with no duplicates (can use union, intersection, difference)

use {} to define sets

```
brca1={"DNA repair", "zinc ion binding", "DNA binding", "ubiqutin"}

# even if you specify a set with repeated elements, python only keeps the unique ones
brca1={"DNA repair", "DNA repair", "zinc ion binding", "DNA binding", "ubiqutin"}

brca2={"DNA repair", "protein binding", "H4 acetylation"}

# union:  |
brca1 | brca2

# intersection:  &
brca1 & brca2

# difference -
brca1 - brca2
brca2 - brca1
```

## dictionaries 

key-value pairs

keys can be string or number (immutable), values can be any type

```
TF_motif = {"SP1" : "GGGCGG",
            "C/EBP" : "ATTGTAG",
            "TBP" : "AGTGATGC"}
TF_motif["SP1"]
print ("SP1 recognizes %s" % TF_motif["SP1"] )
```
check if a key is present:
```
"SP1" in TF_motif
"foo" in TF_motif
```

Modifying dictionaries:
```
# add a new value
TF_motif["AP-1"] = "AACTGTAGT"
# add several new values
TF_motif.update ( {"new1":"AACTGTAGT", 
"new2":"AGTGTAGCT","new3":"GGGGTAGC"} )
# if one of those is an existing key, it updates that one and adds the new ones
TF_motif.update ( {"new4":"AACTGTAGT", "new2":"TTTTT"} )
# modify existing 
TF_motif["AP-1"] = "AACTGTAGTAAA"
# delete a key:
del TF_motif["AP-1"] 
```

Some functions for dictionaries:
```
# length
len TF_motif

# get keys
TF_motif.keys()
# get keys and put them in a list
list(TF_motif.keys())

# get values (in arbitrary order)
TF_motif.values()
list(TF_motif.values())
```

## Ifs and Loops
COLON character is important  
INDENTATION is important  - should be same for each line within the IF loop, and if an ELSE appears, it should be lined up with the IF:
```
# example with n
dna = "atgtangctggatcgtgantcg"
# example without n
dna = "atgtagctggatcgtgatcg"

if 'n' in dna :
  nbases=dna.count('n')
  print("dna seq has %d undefined bases" % nbases)
else :
  print("dna seq has NO undefined bases")
```

Comparison operators (all work for strings and for numbers):
==
!=
<
>
<=
>=

Membership operators: `in` / `not in`

Identity operators: `is` / `is not`.  Similar to, but different from == (and !=). Tests whether objects are in the same location in memory as each other.

Example: using the slice operator ([:]) to create newalphabet means we make a whole new object, so `is` tells us they're not the same

```
alphabet=['a','c','g','t']
newalphabet=alphabet[:]

alphabet==newalphabet 
    # True
alphabet is newalphabet 
    # False
```
but without the slice operator, we didn't create a new object:
```
newalphabet2=alphabet
alphabet is newalphabet2 
    # True

# yikes! this append operation changes alphabet AND newalphabet2 (but not newalphabet), because newalphabet2 is not an independent copy of alphabet:
alphabet.append('n')
alphabet is newalphabet2 
```

More than two alternatives - use `elif`:

```
# example with n
dna = "atgtangctggatcgtgantcg"
# example without n
dna = "atgtagctggatcgtgatcg"
# example without n or c
dna = "atgtagtggatgtgatg"

if 'n' in dna :
  nbases=dna.count('n')
  print("dna seq has %d undefined bases" % nbases)
elif 'c' in dna :
  cbases=dna.count('c')
  print("dna seq has %d C bases" % cbases)
else :
  print("dna seq has NO undefined bases and no Cs")
```

Logical operators:  
`and`  `or`  `not`

```
# example with n
dna = "atgtangctggatcgtgantcg"
# example with n and N
dna = "atgtangctggatcgtgantcgNNNNN"
# example without n
dna = "atgtagctggatcgtgatcg"
# example without n or c
dna = "atgtagtggatgtgatg"

if 'n' in dna or 'N' in dna :
  nbases=dna.count('n') + dna.count('N')
  print("dna seq has %d undefined bases" % nbases)
```

Loops:  `while` loops, `for` loops

While loop (again, indentation is important):
```
dna = "atgtagctgtgatcgtgatcgt"

# if find fails to find anything, it returns -1
pos=dna.find("gt")
while pos>-1 :
  print("found GT at position %d" % (pos+1) )
  pos=dna.find("gt", pos+1)
```

For loop:
```
motifs=["attccgt", "agggggtttttcg", "gtagc"]
for m in motifs:
  print (m,len(m))
```

Use range for numbers:
```
# range(4 returns 0,1,2,3)
for i in range(4):
  print(i)

# this starts at 1, counts by 2, ends when i is still <10 (9 in this case)
for i in range(1,10,2):
  print(i)
```

Check for invalid amino acids:
```
protein="SDVGANNNXCHXP"

## using a list for acceptable_aas - either one of these solutions works:
acceptable_aas =["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

for i in range(len(protein)):
  aa = protein[i]
  if not aa in acceptable_aas:
    print("bad aa %s at position %d" % (protein[i], (i+1)))
# also fine: 
# if aa not in acceptable_aas:

## or, can using a character for acceptable_aas
acceptable_aas_2 = "ARNDCQEGHILKMFPSTWYV"
for i in range(len(protein)):
  aa = protein[i]
  if not aa in acceptable_aas_2:
    print("bad aa %s at position %d" % (protein[i], (i+1)))
# also fine: 
# if aa not in acceptable_aas_2:

## use 'break' to exit the loop at the very first invalid character:
for i in range(len(protein)):
  aa = protein[i]
  if not aa in acceptable_aas:
    print("invalid protein!")
    break

## 'continue' moves to the next iteration of the loop, rather than breaking out of the entire loop (like 'next' in perl)
protein="SDVGANNNXCHXP"
corrected_protein=""
for i in range(len(protein)):
  aa = protein[i]
  if not aa in acceptable_aas:
    continue
  corrected_protein=corrected_protein+aa

print("          protein: %s" % protein)
print("corrected_protein: %s" % corrected_protein)

## else clauses that go with for or while loops. 
# use for a statement you want to run after a loop, except when loop ends with a BREAK statement
# example - find all prime numbers smaller than a given integer (N)
# the else is useful because I want to pay attention to times I went through the whole loop without satisfying the condition
N=10
for y in range(2,N):
  for x in range(2,y):
    if y % x == 0:
      # not a prime, if it's evenly divisible by something
      print("the number %s is divisible by %s" % (y,x))
      break
  else:
    print("the number %d is a prime" % y)
```

`pass` is a placeholder that does nothing. sometimes helps with readability, or as a temporary placeholder when you're developing code
```
if motif not in dna:
  pass
else:
  do_some_function(motif, dna)
```

## Functions

general syntax:
```
def function_name(input arguments) :
    "string documenting the function" # (optional, but nice, as you see it when you do "help(function_name)" )
    function_code_block
    return output
```
again, idents are very important and define where the function starts and ends

```
def gc(dna) :
    """
    compute GC percentage
    more info on that
    """ 
    nbases=dna.count("n") + dna.count("N")
    gcpercent=float( dna.count("c") + dna.count("C") + dna.count("g") + dna.count("G") )*100.0 / (len(dna)-nbases)   
    return gcpercent

dna1 = "atgtagctgcgatcgtgatcg"
gc(dna1)
help(gc)
```

Scope:  global versus local

Types of function:  
Boolean : returns True/False

Task : check for stop codon
```
def check_for_stops(dna) :
    "a function to check a DNA sequence for in-frame stop codons"
    any_stop_codon = ["taa","tga","taa"]
    has_stop = False
    for startPos in range(0,len(dna),3) :
        thisCodon = dna[ startPos : (startPos+3) ].lower()
        if thisCodon in any_stop_codon : 
            has_stop = True
            break
    return(has_stop)

dna1 = "atgtagctgcgatcgtgatcg"
dna2 = "atgtagctgcgatcgtcatcg"
check_for_stops(dna1)
check_for_stops(dna2)
```




Task : check for stop codon, but allow specification of reading frame 0,1,2 (with default 0)
```
def check_for_stops(dna, frame=0) :
    "a function to check a DNA sequence for in-frame stop codons"
    any_stop_codon = ["taa","tga","taa"]
    has_stop = False
    for startPos in range(frame,len(dna),3) :
        thisCodon = dna[ startPos : (startPos+3) ].lower()
        if thisCodon in any_stop_codon : 
            has_stop = True
            break
    return(has_stop)

dna1 = "atgtagctgcgatcgtgatcg"
check_for_stops(dna1,0)
check_for_stops(dna1,1)
check_for_stops(dna1,2)
```

above we passed the argument by position. we can also do it by name:
```
check_for_stops(frame=2,dna=dna1)
```

reverse-complement:

```
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

revcomp(dna1)
```

list comprehensions: a concise way to create a list. A bit like lapply in R.  General syntax:
```
new_list = [operation(i) for i in old_list if filter(i)]
```

`split` and `join` string methods:
```
mysentence = "a few words in a sentence"
mysentence.split()
mysentence.split("in")

mylist = ["a","b","c","d"]
"-".join(mylist)
"".join(mylist)
```
join to me is unintuitive - I'd think it would be a list method not a string method


functions with variable number of arguments - use a * in the definition, and it'll capture remaining arguments as a list
```
def newFunc(arg1, arg2, *theRest) :
    print("arg1 %s" % arg1)
    print("arg2 %s" % arg2)
    print("the rest", theRest)
    return

newFunc("hello","goodbye")
newFunc("hello","goodbye","more1")
newFunc("hello","goodbye","more1","more2")
```

## Modules and packages

A `module` - some functions collected in a file. Has extension `.py`.  Example: dnautil.py

When we want to use it, we do something like this (without the `.py` extension):
```
import dnautil
```

where will python look for modules?
- built-ins
- current working dir
- directory where python is installed
- in a path - a colon-separated list of file paths, stored in PYTHONPATH

Or, we can specify location directly within the script

Show paths, using built-in sys module:
```
import sys
sys.path
```
Use sys.path to extende path within a script:
```
# this doesn't work at first if I'm not in the right dir
import dnautil

# add to path:
import sys
sys.path.append("/Volumes/malik_h/user/jayoung/git_more_repos/MalikLab_bioinformaticsResources/janets_NOTES_forMyself/python")

# now it does work:
import dnautil

# to use functions from modules, we have to use the module name:
dna1 = "atgtagctgcgatcgtgatcg"
dnautil.gc(dna1)
```
alternatively, we load the module like this, and we don't need to specify the module name for each function:
```
import sys
sys.path.append("/Volumes/malik_h/user/jayoung/git_more_repos/MalikLab_bioinformaticsResources/janets_NOTES_forMyself/python")

from dnautil import *
dna1 = "atgtagctgcgatcgtgatcg"
gc(dna1)

# we could have only imported selected functions:
from dnautil import gc,revcomp
```

Sometimes you have >1 package that contains functions with the same name, so you want to use the `dnautil.gc()` way of calling the function, rather than importing all functions in a package (using `from dnautil import *`) and calling them with the shorter name

Packages   
- a group of modules.
- Example:  module name `A.B` means it is submodule `B` in package `A`  
- packages are in DIRECTORIES that contain a special file called `__init__.py`. Can be empty (usually is), but it tells us that the directory is a python package, allowing us to import it.

Example structure for a package dir:
```
bioseq/
  __init__.py
  dnautil.py
  rnautil.py
  proteinutil.py
```

Example structure for a package dir, with sub-packages:
```
bioseq/
  __init__.py
  dnautil.py
  rnautil.py
  proteinutil.py
  fasta/
    __init__.py
    fastautil.py
  fastq/
    __init__.py
    fastqutil.py
```

Could then do 
```
# this
import bioseq.dnautil
bioseq.dnautil.gc("dna")

# or this:
from bioseq import dnautil
dnautil.gc("dna")

# or this to import only one function:
from bioseq.fasta.fastautil import fastaseqread
fastaseqread("fasta")
```

## File I/O

`open` function is for reading and writing

```
# read mode (it's actually the default)
f=open("myfile","r")

# write
f=open("myfile","w")

# append
f=open("myfile","a")

# test for existence
try:
  f = open("myfile","r")
except IOError:
  print ("the file doesn't exist")
```

`IOError` is a special variable returned by `open`

```
f=open("myfile","r")
for line in f:
    print(line)

# alternative:
# if we do this right after the for loop, there's nothing left to read, because file handle is at the end of the file
f.read()

# to re-start a file:
# f.seek(offset, from_what)
f.seek(0)  # resets position to the start

f.read() # slurps the entire file in

f.readline() # reads the first line

f.write(string)   # writes a line to a file (must be in append/write). also returns the number of characers written

# close the file:
f.close()
```

read fasta file and make a dictionary (ignore descriptions):
```
# open file
mySeqFile = "janets_NOTES_forMyself/python/test.fa"

try:
    f = open(mySeqFile, "r")
except IOError: 
    print("warning: file %s does not exist" % mySeqFile)

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

# show all the things in the dictionary
for eachID in mydict.keys():
    print("seqID",eachID,"sequence",mydict[eachID])

# or, fancier:
for name,seq in mydict.items():
    print(name,seq)
```

## command-line arguments

```
python3 processfasta.py myfile.fa
```

inside python we use:
```
import sys
print(sys.argv)
  # showing it's a list with the two items we specified on the command-line
```

see `processfasta.py`


Using `getopt`

Maybe we want an option called -l to specify length threshold

`sys.exit()`  exits

can separate two commands with semicolon

can use `stdin` stream

output (unless otherwise specified) goes to `stdout` (capture using `1>`)

there's also `stderr` (capture using `2>`)

The `sys` module can help get stdin

```
sys.stdin.read()
    # this will wait for input.  
    # user presses ctrl-D to end the input

sys.stdout.write("hello")
    # will also return the length of the string

sys.stderr.write("warning")
    # will also return the length of the string
```


calling other programs:  `call/execute` from the `subprocess` module

```
import subprocess
# call expects a list as input:
subprocess.call(["ls", "-l"])
    # returns 0 as well as the output
subprocess.call(["tophat", "genomeIndex", "forwardReads.fq", "reverseReads.fq"])
```

## [Biopython](https://biopython.org)

Started in 1999

Parsing, online database access, interfaces to various common programs.

```
import Bio
print(Bio.__version__)
```

See `testBiopython.py`

[Tutorial/cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html)

[FAQs](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec6)


## More functions

replace for strings:
```
mystring="hello"
mystring.replace("l","")
```
re.sub : regular expression replace
```
import re
s = "123123"
s = re.sub('23$', 'penguins', s)
print(s)
```

## Installations

My mac laptop did not have pip on it, so I do this;

First I install the latest version of python, using `brew`:


This installs a SECOND copy of python3 in a non-default location, alpong with pip and setuptools
```
brew install python
```
Messages include:
```
==> python@3.11
Python has been installed as
  /usr/local/bin/python3

Unversioned symlinks `python`, `python-config`, `pip` etc. pointing to
`python3`, `python3-config`, `pip3` etc., respectively, have been installed into
  /usr/local/opt/python@3.11/libexec/bin

You can install Python packages with
  pip3 install <package>
They will install into the site-package directory
  /usr/local/lib/python3.11/site-packages
```

I did this (not sure if it actually helped) to make symlinks again:
```
brew unlink python && brew link python
```

I want to add `/usr/local/opt/python@3.11/libexec/bin` to my PATH so that I can use `pip` to mean `pip3` and `python` to mean `python3` (in each case, this new version I just installed). So I modify `.profile` to add the following line:
```
export PATH="/usr/local/opt/python@3.11/libexec/bin:$PATH"
```

Now pip and python are both in my PATH
```
which python
    # /usr/local/opt/python@3.11/libexec/bin/python
which pip
    # /usr/local/opt/python@3.11/libexec/bin/pip
```

```
which python3
    # /usr/local/bin/python3
ls -l /usr/local/bin/python3
lrwxr-xr-x  1 jayoung  admin  42 Dec 28 15:31 /usr/local/bin/python3 -> ../Cellar/python@3.11/3.11.6_1/bin/python3
```

I think this is what I'd like the header of my phython scripts to have, on my Mac laptop at least:
`/usr/local/bin/python3`

Now I can install Biopython:
```
pip install biopython
   # Installing collected packages: numpy, biopython
   # Successfully installed biopython-1.82 numpy-1.26.2
```

where did it get installed?
```
# show which packages are installed and where they are:
python3 -m pip list -v
```

on my mac laptop, packages are here: `/usr/local/lib/python3.11/site-packages`

So I want to add that to PYTHONPATH - I add this to my .profile file:
```
#### JY adding to (or creating) PYTHONPATH:
if [ -z ${PYTHONPATH+x} ] 
then
    echo "        making new PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages
else
    echo "        adding to PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages:${PYTHONPATH}
fi
```


