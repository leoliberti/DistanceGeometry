#!/usr/bin/env python
import sys
from dg import *

inputFile = sys.argv[1]
dmdgp = instanceFromFile(inputFile)
G = dmdgp.graph

x = SDP(G)
print LDE(G, x)

