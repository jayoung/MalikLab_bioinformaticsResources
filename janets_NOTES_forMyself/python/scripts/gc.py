#!/usr/local/bin/python3

# dna = "atgtagctggatcgtgatcg"

dna = input("Enter your DNA seq please:")

numC = dna.count("c")
numG = dna.count("g")

gcPercent = 100*(numC+numG) / len(dna)

# print("gcPercent is", gcPercent, "%")

### with formatting
print("gcPercent is %5.1f %%" % gcPercent)

print("done")
