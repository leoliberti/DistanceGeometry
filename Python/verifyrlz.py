#!/usr/bin/env python

## name: verifyrlz.py
## author: Leo Liberti
## purpose: verify MDE/LDE of a realization for a given DGP instance
## source: python 2.7
## history: 160613 derived from isomap.py

import sys
from dg import *

maxNLPitn = 50

####################### MAIN ########################

if len(sys.argv) < 3:
    print "error, need instance file and realization file on cmd line"
    print "instance file formats: .gph, .nmr, .dimacs, .dgsol, .txt, .dat"
    print "realization file format: dgsol's .sol output"
    exit(1)

# set number of dimensions
NumberOfDimensions = 3

# read graph from data file
inputFile = sys.argv[1]
rlzFile = sys.argv[2]
G = Graph(edgesFromFile(inputFile))
x = realizationFromFile(G, rlzFile)
filename = os.path.basename(inputFile)

print filename, G.V.card, len(G.E), MDE(G,x), LDE(G,x)

#print x
#G.plot(x)
