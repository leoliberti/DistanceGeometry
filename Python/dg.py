#####################################
# filename: dg.py
# author: Leo Liberti
# purpose: module for doing distance geometry
# source: Python 2.7
# history: 160524 union of dg_isomap.py, formulations.py, bp.py
# CONVENTION: realizations are either lists of lists, or matrices;
#             in either case, x[v,k] is the k-the component of vertex v
# requirements: Python >= 2.7.11,
#               pip install numpy scipy matplotlib algopy cvxopt picos
#               pyipopt, mosek, cplex (python $i install)
###################################

############# SECT:IMPORTS #############
# basic
import sys         # system calls (.exit, .argv)
import re          # regular expression library, used in reading input files
import os.path
import time
import math        # standard math functions 
import numpy as np # the numeric library, with linear algebra calls
#from numpy import *
import scipy
from scipy import optimize # to use the optimization solvers in scipy 
from collections import OrderedDict # used to nice-print adjacency lists
import heapq       # priority queue for non-recursive BP? [faster]
import Queue as Q  # priority queue for non-recursive BP? [thread-safe]
#from heapq_showtree import show_tree
import itertools   # find initial clique in DVOP order greedy alg
import operator    # operator.itemgetter() to find max in tuple components
import types
import warnings
#import nlopt
#import ad
import copy
from functools import total_ordering # used in BP

# plots
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

# CPLEX
import cplex
from cplex.exceptions import CplexError

# local NLP
#import pyipopt       # IPOpt API
import functools     # automatic differentiation - Jacobian / Hessian (helper)
import algopy        # automatic differentiation - Jacobian / Hessian 
#import dgopt         # the DGSol library ported to Python through F2PY

# SDP
import cvxopt as cvx # the CvxOpt library, with various solvers including SDP
import picos as pic  # the modelling layer for CvxOpt

# options
np.set_printoptions(precision=2, suppress=True) # print precision for numpy

############# SECT:GLOBALS ###############

# set number of dimensions
global NumberOfDimensions
NumberOfDimensions = 3

# this is a technical hack to plot in 3D and update the plot after mouse ops
#   no need to change / update this; see
#   stackoverflow.com/questions/10374930/matplotlib-annotating-a-3d-scatter-plot
global labels_and_points
labels_and_points = []

# this is a tolerance for zero, used in next* computations
global myZero
myZero = 1e-6

# this is an upper bound on weights
global myInf
myInf = 1e+25

# iterations of DDP iterative method
global DDPItn
DDPItn = 5

# maximum iteration for local NLP solvers (<0 to accept solver's default)
nlpMaxIter = -1

# SDP solver tolerance (<0 to accept solver's default)
solverTolerance = 1.0e-3

# number of samplings in Barvinok's naive algorithm (ISCO16=5)
BarvinokIterations = 500

############# SECT:CLASSES #############

## this class defines a vertex
##   a vertex has name, rank and an optional Euclidean position
class Vertex:
    # declare and initialize member attributes (rk start from 0)
    def __init__(self, nm, rk=-1):
        # a vertex is name and rank, and other optional fields
        self.name = nm
        self.rank = rk
        self.pos = list()
    # printing vertices
    def __str__(self):
        return str(self.name)


## this class defines an edge or arc
##   an edge has a tail, a head, a lower/upper weights, a directed flag
##   optionally, also some chemical information for molecule graphs
##   for non-interval weights, wL represents the weight and wU is ignored
class Edge:
    # declare and initialize member attributes
    def __init__(self, v1, v2, w1, w2, givenDir):
        self.tail = v1
        self.head = v2
        self.wL = w1   
        self.wU = w2
        self.interval = 0
        if abs(w2 - w1) > myZero:
            self.interval = 1
        self.directed = givenDir
        # optional fields, only useful for .nmr files
        self.tlatom = ""
        self.hdatom = ""
        self.amino1 = ""
        self.amino2 = ""
    # printing edges
    def __str__(self):
        if self.directed == True:
            link = "->"
        else:
            link = "--"
        return str(self.tail) + link + str(self.head)
    # initialize the additional data structures in NMR files
    def nmr(self, atom1, atom2, am1, am2):
        self.tlatom = atom1
        self.hdatom = atom2
        self.amino1 = am1
        self.amino2 = am2

## this class defines a set of vertices
##   vertices are stored by name or by rank; this class provides the mappings
class Vertices:
    # initialize vertices from a set of edges
    def __init__(self, E):
        # create a set of vertex names
        self.names = set()
        # and fill it from the edges
        for e in E:
            self.names.add(e.tail.name)
            self.names.add(e.head.name)
        self.names = sorted(self.names, key=int)
        # create maps to find vertices from rank and name
        self.byRank = dict()
        self.byName = dict()
        i = 0   # this marks the fact that vertex ranks start from 0
        maxrank = 0
        # fill these maps
        for vn in self.names:
            v = Vertex(vn)
            if v.rank == -1:
                v.rank = i
            if v.rank > maxrank:
                maxrank = v.rank
            if v.rank in self.byRank:
                print "dg.py:Vertices: two vertices have same (pre-assigned) rank"
                exit('abort')
            if v.rank in self.byName:
                print "dg.py:Vertices: two vertices have same name"
                exit('abort')
            self.byRank[v.rank] = v
            self.byName[vn] = v
            i = i + 1
        # set the correct vertex rank within the edge set
        for e in E:
            e.head.rank = self.byName[e.head.name].rank
            e.tail.rank = self.byName[e.tail.name].rank
        # cardinality of the vertex set
        self.card = i
        if self.card != maxrank + 1:
            print "dg.py:Vertices: vtxrk inconsistency: |V| =", self.card, "!=", maxrank+1, "= maxrk+1"
            exit('abort')
    # print vertex sets
    def __str__(self):
        return " ".join(self.names)

## we don't define a set of edges, since it's simply a list of edges
##   it's the first data structure we read from input files
        
## this class defines a weighted graph
##   a graph has a Vertices object and a list of edges
##   it computes adjacency lists, which are somewhat cumbersome to use
##   example:
##     for u in range(n):
##       for v in G[G.V.byRank[u].name]:
##         # do something with (edge) G[G.V.byName[u].name][v]
class Graph:
    # declare and initialize member attributes
    def __init__(self, E):
        # vertices
        self.V = Vertices(E)
        # edges
        self.E = E
        # adjacency lists (implemented as a dict:V->(dicts:V->edge))
        #   maps v to dictionary mapping u to directed edge (v,u) and,
        #   if graph undirected, also store inverse edge in head list
        self.adj = dict()
        for v in self.V.names:
            self.adj[v] = dict()
        for e in E:
            self.adj[e.tail.name][e.head.name] = e
            if e.directed == 0:
                # inverse edge appearing in head list is directed as -1
                einv = Edge(e.head, e.tail, e.wL, e.wU, -1)
                einv.nmr(e.tlatom, e.hdatom, e.amino1, e.amino2)
                self.adj[e.head.name][e.tail.name] = einv
    # access the v-the adjacency list as G[v]
    def __getitem__(self, vname):
        return self.adj[vname]
    # return {u,v} if it is an edge, None otherwise
    def getEdge(self, u,v):
        try:
            e = self.adj[self.V.byRank[u].name].get(self.V.byRank[v].name)
        except KeyError, errmsg:
            e = None
        return e
    def getEdgeByName(self, u,v):
        return self.adj[u].get(v)
    def isEdge(self, u,v):
        try:
            e = self.adj[self.V.byRank[u].name].get(self.V.byRank[v].name)
        except KeyError, errmsg:
            e = None
        if (e != None):
            return True
        else:
            return False
    def isEdgeByName(self, u,v):
        e = self.adj[u].get(v)
        if (e != None):
            return True
        else:
            return False
    # return True iff {u,v} is an edge, where u,v are ranks
    def getEdgeWeight(self, u,v):
        try: 
            e = self.adj[self.V.byRank[u].name].get(self.V.byRank[v].name)
        except KeyError, errmsg:
            e = None
        if (e != None):
            return (e.wL, e.wU)
        else:
            return e
    # return the weight interval of the edge {u,v} where u,v are ranks
    # print the adjacency lists
    def __str__(self):
        adj = OrderedDict(sorted(self.adj.items(), key = lambda t : int(t[0])))
        ret = str()
        for v in adj:
            ret = ret + v + ":"
            for u in adj[v].values():
                s = " {0:s}({1:.2f})".format(u.head.name, round(u.wL, 2))
                ret = ret + s
            ret = ret + "\n";
        return ret
    # return a sparse version of the weighted adjacency matrix of the graph
    def sparseWeightedAdjMatrix(self):
        I = list()
        J = list()
        W = list()
        for e in self.E:
            i = e.tail.rank
            j = e.head.rank
            w = e.wL
            W.append(w)
            I.append(i)
            J.append(j)
            if e.directed == 0:
                W.append(w)
                I.append(j)
                J.append(i)
        return scipy.sparse.csr_matrix(scipy.sparse.coo_matrix((W,(I,J)), shape = (self.V.card,self.V.card)))
    # return the weighted adjacency matrix of the graph
    def weightedAdjMatrix(self):
        n = self.V.card
        A = myInf * np.ones((n,n))
        for i in range(n):
            A[i,i] = 0
        for e in self.E:
            A[e.tail.rank,e.head.rank] = e.wL
            A[e.head.rank,e.tail.rank] = e.wL
        return A
    # turn all wL=wU to (wL-err,wU+err) where err=percent*|wU|,
    #   return number of changed edges
    def intervalWeights(self, percent = 0.1):
        changed = 0
        for e in self.E:
            if abs(e.wU - e.wL) < myZero:
                err = percent*abs(e.wU)
                e.wL = e.wL - err
                e.wU = e.wU + err
                e.interval = 1
                changed = changed + 1
        return changed
    # plot the graph in 1D, 2D or 3D using the realization x
    #   R is an instance vertex rank order
    #   if R = [], no edge is drawn
    #   if len(R) == 1, the whole edge set is drawn
    #   if len(R) >= G.V.card, only the edges between R[i] and R[i-1] are drawn
    def plot(self, x, R = []):
        K = NumberOfDimensions
        n = self.V.card
        if type(x) is dict:
            x = dict2array(x)
        fig = plt.figure()
        if K < 1:
            K = 1
        elif K > 3:
            K = 3
            x = x[:,0:K]
        if K == 1:
            # one-dimensional plots on a line y=1 in the xy plane
            ax = fig.add_subplot(111)
            # points and labels
            y = np.ones(len(x[:,0]))
            ax.plot(x[0,:], y, c='r', marker='o', linestyle='None')
            for label,lx,ly in zip(self.V.names,x[:,0], y):
                ax.annotate(label, xy=(lx,ly), xytext=(7,-3), textcoords='offset points')
            # edges as curved segments above and below
            if len(R) > 0:
                # only print edges if R not empty
                i = 0
                for e in self.E:
                    tl = e.tail.rank
                    hd = e.head.rank
                    if len(R) >= n:
                        # a (re)order is given, check e is in the order
                        edgeinorder = False
                        try:
                            tlidx = R.index(tl)
                            if tlidx == 0:
                                if R[1] == hd:
                                    edgeinorder = True
                            elif tlidx == n-1:
                                if R[n-2] == hd:
                                    edgeinorder = True
                            elif R[tlidx-1] == hd or R[tlidx+1] == hd:
                                edgeinorder = True
                        except ValueError:
                            pass
                        if not edgeinorder:
                            # e not in R, don't print it
                            continue
                    updown = 1.5
                    if i % 2 == 1:
                        updown = 0.5
                    verts = [(x[tl,0], 1.), (0.5*(x[tl,0]+x[hd,0]), updown), (x[hd,0], 1.)]
                    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
                    path = Path(verts, codes)
                    patch = patches.PathPatch(path, facecolor='none', edgecolor='blue', lw=1)
                    ax.add_patch(patch)
                    i = i + 1
            # axes
            epsx = 0.1*(max(x[:,0])-min(x[:,0]))
            xmin = min(x[:,0]) - epsx
            xmax = max(x[:,0]) + epsx
            ax.axis([xmin, xmax, 0, 2])
        elif K == 2:
            # two-dimensional plots on the xy plane
            ax = fig.add_subplot(111)
            # points and labels
            ax.plot(x[:,0], x[:,1], c='r', marker='o', linestyle='None')
            for label,lx,ly in zip(self.V.names,x[:,0],x[:,1]):
                ax.annotate(label, xy=(lx,ly), xytext=(-6,6), textcoords='offset points')
            # edges as straight segments
            if len(R) > 0:
                # only print edges if R not empty
                for e in self.E:
                    tl = e.tail.rank
                    hd = e.head.rank
                    if len(R) >= n:
                        # a (re)order is given, check e is in the order
                        edgeinorder = False
                        try:
                            tlidx = R.index(tl)
                            if tlidx == 0:
                                if R[1] == hd:
                                    edgeinorder = True
                            elif tlidx == n-1:
                                if R[n-2] == hd:
                                    edgeinorder = True
                            elif R[tlidx-1] == hd or R[tlidx+1] == hd:
                                edgeinorder = True
                        except ValueError:
                            pass
                        if not edgeinorder:
                            # e not in R, don't print it
                            continue
                    ax.plot([x[tl,0], x[hd,0]], [x[tl.rank,1], x[hd,1]], c='b')
            # axes
            epsx = 0.1*(max(x[:,0])-min(x[:,0]))
            epsy = 0.1*(max(x[:,1])-min(x[:,1]))
            xmin = min(x[:,0]) - epsx
            xmax = max(x[:,0]) + epsx
            ymin = min(x[:,1]) - epsx
            ymax = max(x[:,1]) + epsx
            ax.axis([xmin, xmax, ymin, ymax])
        elif K == 3:
            # three-dimensional plots in xyz space
            ax = fig.add_subplot(111, projection='3d')
            # points and labels
            ax.scatter(x[:,0], x[:,1], x[:,2], c='r', marker='o')
            for label,lx,ly,lz in zip(self.V.names,x[:,0],x[:,1],x[:,2]):
                x2, y2, _ = proj3d.proj_transform(lx,ly,lz, ax.get_proj())
                lb = ax.annotate(label, xy=(x2,y2), xytext=(-6,6), textcoords='offset points')
                labels_and_points.append((lb, lx, ly, lz))
            # edges as straight segments
            if len(R) > 0:
                for e in self.E:
                    tl = e.tail.rank
                    hd = e.head.rank
                    if len(R) >= n:
                        # a (re)order is given, check e is in the order
                        edgeinorder = False
                        try:
                            tlidx = R.index(tl)
                            if tlidx == 0:
                                if R[1] == hd:
                                    edgeinorder = True
                            elif tlidx == n-1:
                                if R[n-2] == hd:
                                    edgeinorder = True
                            elif R[tlidx-1] == hd or R[tlidx+1] == hd:
                                edgeinorder = True
                        except ValueError:
                            pass
                        if not edgeinorder:
                            # e not in R, don't print it
                            continue
                    ax.plot([x[tl,0], x[hd,0]], [x[tl,1], x[hd,1]], [x[tl,2], x[hd,2]], c='b')
            # axes
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            def update_position(e):
                for lb, lx, ly, lz in labels_and_points:
                    x2, y2, _ = proj3d.proj_transform(lx, ly, lz, ax.get_proj())
                    lb.xy = x2,y2
                    lb.update_positions(fig.canvas.renderer)
                fig.canvas.draw()
            fig.canvas.mpl_connect('motion_notify_event', update_position) 
        plt.show()

## this class defines a DDGP instance
##   the order R (on vertex names, not ranks!) is either given or left empty
##   has to be given in DMDGP cases, otherwise natural order assumed DMDGP
##   DDGP cases: if not given it's computed in FPT time (2^K) from DDGPorder()
class DDGP:
    def __init__(self, G, R, ddgptype = "none"):
        self.ordertype = "order"
        self.pred = {}
        K = NumberOfDimensions
        if G.V.card <= K:
            print "dg.py:DDGP(): only", G.V.card, "vertices in", K, "dimensions, need more"
            exit('abort')
        # ddgptype can be "none" (no vtx ords specified), "ddgp", "dmdgp"
        #   ddgp = predecessors may not be contiguous, R is total order on G.V
        #   dmdgp = contiguous predecessors, R is an order or re-order on G.V
        if len(R) == 0:
            if ddgptype == "ddgp":
                # find a DDGP order
                (R,alpha,flag,rd,rdinv) = DDGPOrder(G, K)
                if flag == False:
                    exit('dg.py:DDGP(): instance is not DDGP')
                else:
                    print "dg.py:DDGP(): found DDGP order", R
            elif ddgptype == "dmdgp":
                # if R is not given, assume the "natural order" on V
                R = [G.V.byRank[v].name for v in range(G.V.card)]
        if ddgptype != "none":
            # initialize predecessor dict: vrank->list_of_predranks
            self.pred = {G.V.byName[r].rank:[] for r in R}
            # first K+1 vertices
            for i in range(1,K+1):
                for j in range(i):
                    self.pred[G.V.byName[R[i]].rank].append(G.V.byName[R[j]].rank)
        if ddgptype == "ddgp":
            # in DDGP, R is NOT a re-order
            # r is a vertex rank, but its rank in the order R is i
            for i,r in enumerate(R):
                # the following is NOT equivalent to i,r in enumerate(R[K+1:])
                #   since i is supposed to start from K+1
                if i >= K+1:
                    count = 0
                    for v in G[G.V.byName[r].name]:
                        s = G.V.byName[v].name
                        sinR = True
                        try:
                            j = R.index(s) # rank of the vtxrank s in order R
                        except:
                            # we should only worry about vertices in R
                            sinR = False
                        if sinR:
                            if j < i:
                                self.pred[G.V.byName[r].rank].append(G.V.byName[s].rank)
                                count = count + 1
                    if count < K:
                        print "dg.py:DDGP(): not DDGP, vtx", r, "has only", count, "<", K, "predecessors"
                        exit('abort')
        elif ddgptype == "dmdgp":
            for i,r in enumerate(R):
                if i >= K+1:
                    for j in range(K):
                        if R[i-j-1] == r:
                            self.ordertype = "reorder"
                            # don't bother with a predecessor order, do it
                            # later in BP
                            self.pred = {}
                            break
                        else:
                            if self.ordertype == "order":
                                self.pred[G.V.byName[r].rank].append(G.V.byName[R[i-j-1]].rank)
        self.order = R
        self.graph = G
        self.type = ddgptype

    def ordCheck(self):
        # check order and instance type consistency
        K = NumberOfDimensions
        ddgptype = self.type
        if len(self.order) > 0:
            # some order is given
            dmdgpord = isDMDGPInstance(self, K)
            ret = True
            if ddgptype == "ddgp" and dmdgpord == True:
                print "dg.py:DDGP::ordCheck: DDGP instance with DMDGP order, going on"
            if ddgptype == "dmdgp" and dmdgpord == False:
                print "dg.py:DDGP::ordCheck: DMDGP instance but order is not DMDGP"
                print "   trying to find a DDGP order..."
                (R,alpha,flag,rd,rdinv) = DDGPOrder(self.graph, K)
                if flag == False:
                    exit("dg.py:DDGP::ordCheck: instance not even DDGP, abort")
                else:
                    print "dg.py:DDGP::ordCheck: relabelling instance as DDGP"
                    ret = False
                    self.__init__(self.graph, R, "ddgp")
        else:
            # no order is given
            ret = False
            self.__init__(self.graph, self.order, "ddgp")
        return ret

    # write the instance in some format frm
    # if relabel = True relabel the vertices according to self.order
    def write(self, outputFile, frm = "txt", relabel = False):
        R = self.order
        if len(R) == 0 and relabel == True:
            print "dg.py:DDGP::write(): no relabelling without an order"
            relabel = False
        G = self.graph
        n = G.V.card
        if len(R) > n and relabel == True:
            print "dg.py:DDGP::write(): no relabelling using a reorder"
            relabel = False
        if len(R) == 0 and relabel == False:
            R = [G.V.byRank[i].name for i in range(n)]
        if relabel == True:
            # relabel according to given order
            lbl = {nm:i for i,nm in enumerate(R)}
        else:
            # identity
            lbl = {nm:nm for nm in R}
        f = open(outputFile, "w")
        if frm == "txt":
            print >> f, "# dg.py: DGP .txt-formatted instance", outputFile
            for e in G.E:
                print >> f, lbl[e.tail.name], lbl[e.head.name], e.wL, e.wU, e.interval
        elif frm == "nmr":
            #print >> f, "# dg.py: DGP .nmr-formatted instance", outputFile
            for e in G.E:
                print >> f, lbl[e.tail.name], lbl[e.head.name], e.wL, e.wU, "N N ALA ALA"
        elif frm == "gph":
            print >> f, "# dg.py: DGP .gph-formatted instance", outputFile
            print >> f, "Undirected;"
            print >> f, "Vertices:"
            for r in R:
                print >> f, lbl[r]
            print >> f, "Edges:"
            for e in G.E:
                print >> f, lbl[e.tail.name], lbl[e.head.name], e.wL, e.wU #nonconformant
                #print >> f, lbl[e.tail.name], lbl[e.head.name], e.wL, "1"#gph-conformant
        else:
            print "dg.py:DDGP::write: format", frm, "unknown"
            exit('abort')
        f.close()

    def subInstance(self, S):
        # return a new instance corresponding to the induced subgraph S
        G = self.graph
        R = [G.V.byRank[i].name for i in S]
        return DDGP(inducedGraph(G, S), R, "none")
        
## this class used to prettyprint list of floats
##   use as: print map(floatPrint, myList)
class floatPrint(float):
    def __repr__(self):
        return "%0.2f" % self

        
############# SECT:INSTANCE_INPUT ############

## read input file
##   this a dispatcher procedure (not used since 170503)
def edgesFromFile(inputFile, frm = ""):
    E = list()
    if frm == "":
        frm = fileFormat(inputFile)
    if frm == "gph":
        (E,R) = readGphFile(inputFile)
    elif frm == "nmr":
        (E,R) = readNMRFile(inputFile)
    elif frm == "dimacs":
        (E,R) = readDIMACSFile(inputFile)
    elif frm == "txt":
        (E,R) = readTxtFile(inputFile)
    elif frm == "dat":
        (E,R) = readDatFile(inputFile)
    else:
        sys.exit("{0:s}: error: can't handle format \"{1:s}\" in input file {2:s}".format(sys.argv[0],frm,inputFile))
    return E

## read data file, decide input format
def fileFormat(filename):
    theFormat = "unknown"
    with open(filename) as f:
        for line in f:
            # look for first non-comment line
            line = line.strip()
            if len(line) > 0:
                if line[0] != '#':
                    # recognize format
                    if line == "Undirected;" or line == "Directed;":
                        # the .gph format (see makegraph.c)
                        theFormat = "gph"
                    elif line[0] == 'c' and line[1] == ' ':
                        # the well-known DIMACS format
                        theFormat = "dimacs"
                    elif line[0:5] == 'param':
                        # AMPL .dat file (automatically translated from gph)
                        theFormat = "dat"
                    elif line.lower() == "begin vertices":
                        # Andrea Cassioli's .dmdgp format
                        theFormat = "dmdgp"
                    else:
                        cols = line.split()
                        if len(cols) < 3:
                            theFormat = "unknown"
                        else:
                            m1 = re.search("[0-9]*", cols[0])
                            vtx1 = m1.group(0)
                            m2 = re.search("[0-9]*", cols[1])
                            vtx2 = m2.group(0)
                            m3 = re.search("[0-9]*\.?[0-9]*", cols[2])
                            ew1 = m3.group(0)
                            if len(cols) >= 4:
                                m4 = re.search("[0-9]*\.?[0-9]*", cols[3])
                                ew2 = m4.group(0)
                                if vtx1 == cols[0] and vtx2 == cols[1] and ew1 == cols[2] and ew2 == cols[3]:
                                    # the .nmr format (PDB)
                                    theFormat = "nmr"
                            elif vtx1 == cols[0] and vtx2 == cols[1] and ew1 == cols[2]:
                                # the simplest .txt format [vtx1 vtx2 weight12]
                                theFormat = "txt"
                            else:
                                theFormat = "unknown"
                    # we have the format, stop scanning
                    break
    if theFormat == "nmr":
        print "dg.py:fileFormat(): input format is txt or nmr"
    else:
        print "dg.py:fileFormat(): input format is", theFormat
    return theFormat

## read graph (& vtx order) from input file
def instanceFromFile(inputFile, dgptype = "none", frm = "none"):
    if frm == "none":
        frm = fileFormat(inputFile)
    E = list()
    R = list()
    if frm == "gph":
        (E,R) = readGphFile(inputFile)
    elif frm == "nmr":
        (E,R) = readNMRFile(inputFile)
    elif frm == "dimacs":
        (E,R) = readDIMACSFile(inputFile)
    elif frm == "txt":
        (E,R) = readTxtFile(inputFile)
    elif frm == "dat":
        (E,R) = readDatFile(inputFile)
    elif frm == "dmdgp":
        (E,R) = readDMDGPFile(inputFile)
    else:
        sys.exit("{0:s}: error: can't handle format \"{1:s}\" in input file {2:s}".format(sys.argv[0],frm,inputFile))
    G = Graph(E)
    if frm == "dmdgp":
        dgp = DDGP(G,R, "dmdgp")
    else:
        dgp = DDGP(G,R, dgptype)
    return dgp

## read the DIMACS graph file format
def readDIMACSFile(filename):
    dirFlag = False
    E = list()
    R = list()
    with open(filename) as f:
        for linecount, line in enumerate(f):
            line = line.strip()
            if len(line) > 0:
                # we're only interested in "arc" lines
                if line[0] == 'a':
                    cols = [c for c in line.split() if not '#' in c]
                    if len(cols) >= 5:
                        e = Edge(Vertex(cols[1]), Vertex(cols[2]), float(cols[3]), float(cols[4]), dirFlag)
                        if e.wL > e.wU:
                            print "readDIMACSFile: WARNING: interval weight[", e.wL, ",", e.wU, "] empty, setting to", e.wL
                            e.wU = e.wL
                    elif len(cols) >= 4:
                        e = Edge(Vertex(cols[1]), Vertex(cols[2]), float(cols[3]), float(cols[3]), dirFlag)
                    else:
                        print "readDIMACSFile: ERROR: arc line", linecount, "has < 4 columns"
                        exit('abort')
                    E.append(e)
    return (E,R)

## read data file written in .dmdgp format (format by Andrea Cassioli)
def readDMDGPFile(filename):
    dirFlag = False
    # we read F, then the BP order, then we permute edges and put them in E
    E = list()
    bporder = False
    R = list()
    section = "none"
    with open(filename) as f:
        for i,line in enumerate(f):
            line = line.strip()
            if len(line) > 0:
                if line[0] != '#':
                    # decide what section of the file we're in
                    if line.lower() == "begin vertices":
                        section = "vertices"
                    elif line.lower() == "end vertices":
                        section = "none"
                    elif line.lower() == "begin edges":
                        section = "edges"
                    elif line.lower() == "end edges":
                        section = "none"
                    elif line.lower() == "begin atom_names":
                        section = "atom_names"
                    elif line.lower() == "end atom_names":
                        section = "none"
                    elif line.lower() == "begin residues":
                        section = "residues"
                    elif line.lower() == "end residues":
                        section = "none"
                    elif line.lower() == "begin dihedral_angles":
                        section = "dihedral_angles"
                    elif line.lower() == "end dihedral_angles":
                        section = "none"
                    elif line.lower() == "begin bp_order":
                        section = "bp_order"
                    elif line.lower() == "end bp_order":
                        section = "none"
                    else:
                        # actions for each section
                        if section == "vertices":
                            # TO DO
                            a = False
                        elif section == "edges":
                            # read the (unpermuted) edge
                            cols = line.split()
                            try:
                                cols = cols[0:cols.index('#')]
                            except:
                                pass
                            if len(cols) >= 5:
                                e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[3]), float(cols[4]), dirFlag)
                                if e.wL > e.wU:
                                    print "readDMDGPFile: WARNING: interval weight[", e.wL, ",", e.wU, "] empty, setting to", e.wL
                                    e.wU = e.wL
                            elif len(cols) >= 4:
                                e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[3]), float(cols[3]), dirFlag)
                            else:
                                print "readDMDGPFile: ERROR: line", linecount, "has < 4 columns"
                                exit('abort')                            
                            E.append(e)
                        elif section == "atom_names":
                            # TO DO
                            a = False
                        elif section == "residues":
                            # TO DO
                            a = False
                        elif section == "dihedral_angles":
                            # TO DO
                            a = False
                        elif section == "bp_order":
                            bporder = True
                            # vertices are stored as labels -- mostly, they
                            #   are integers, but they need not be
                            R.append(line.split()[0])
    return (E,R)

## read data file written in AMPL .dat format
## (after line starting "param : E" there is a list of edges formatted as
##    i j w_{ij} I_{ij} )
def readDatFile(filename):
    dirFlag = False
    E = list()
    R = list()
    edgeflag = False
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:                
                if line[0] != '#':
                    if edgeflag:
                        if line[0] == ';':
                            edgeflag = False
                        else:
                            cols = [c for c in line.split() if not '#' in c]
                            if len(cols) >= 6:
                                e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[4]), float(cols[5]), dirFlag)
                            elif len(cols) >= 4:
                                e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[2]), float(cols[2]), dirFlag)
                                if e.wL > e.wU:
                                    print "readDatFile: WARNING: interval weight[", e.wL, ",", e.wU, "] empty, setting to", e.wL
                                    e.wU = e.wL
                            else:
                                print "readDatFile: ERROR: line", linecount, "has < 4 columns"
                                exit('abort')
                            E.append(e)
                    else:
                        if line[0:9].replace(" ","") == 'param:E':
                            edgeflag = True
    return (E,R)

## read data file written by MakeGraph -- also reads vertex order
def readGphFile(filename):
    dirFlag = False
    vtxFlag = True
    E = list()
    R = list()
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] != '#':
                    if vtxFlag == True:
                        # vertex part
                        if line == "Directed;":
                            dirFlag = True
                        elif line == "Undirected;":
                            pass
                        elif line == "Vertices:":
                            pass
                        elif line == "Edges:":
                            vtxFlag = False
                        else:
                            # read vertex order
                            R.append(line)
                    else:
                        # edge part
                        cols = [c for c in line.split() if not '#' in c]
                        if len(cols) >= 4:
                            if cols[3] == "1" and cols[2] != "1" and float(cols[2]) >= 1:
                                # gph inclusion field, not intdist
                                cols[3] = cols[2]
                            e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[2]), float(cols[3]), dirFlag)
                            E.append(e)
                        else:
                            print "readGphFile: ERROR: line", linecount, "has < 4 columns"
                            exit('abort')
    return (E,R)

## read data file written in .nmr format (PDB / Antonio Mucherino)
def readNMRFile(filename):
    dirFlag = False
    E = list()
    R = list()
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] != '#':
                    cols = [c for c in line.split() if not '#' in c]
                    e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[2]), float(cols[3]), dirFlag)
                    # look for optional columns
                    if len(cols) >= 8:
                        e.nmr(cols[4], cols[5], cols[6], cols[7])
                    elif len(cols) >= 7:
                        e.nmr(cols[4], cols[5], cols[6], "")
                    elif len(cols) >= 6:
                        e.nmr(cols[4], cols[5], "", "")
                    elif len(cols) >= 5:
                        e.nmr(cols[4], "", "", "")
                    elif len(cols) < 4:
                        print "readNMRFile: ERROR: line", linecount, "has < 4 columns"
                        exit('abort')
                    E.append(e)
    return (E,R)

## read data file written in .txt format (each line is "vtx1 vtx2 wL [wU]")
def readTxtFile(filename):
    dirFlag = False
    E = list()
    R = list()
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] != '#':
                    cols = [c for c in line.split() if not '#' in c]
                    if len(cols) >= 4:
                        e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[2]), float(cols[3]), dirFlag)
                    elif len(cols) >= 3:
                        e = Edge(Vertex(cols[0]), Vertex(cols[1]), float(cols[2]), float(cols[2]), dirFlag)
                    else:
                        print "readTxtFile: ERROR: line", linecount, "has < 3 columns"
                        exit('abort')
                    E.append(e)
    return (E,R)

## read realization from file (each row has K columns)
def realizationFromFile(G, inputFile):
    n = G.V.card
    K = NumberOfDimensions
    x = np.zeros((n,K))
    f = open(inputFile, "r")
    irk = 0
    inm = 1
    for line in f:
        v = line.split()
        # assumption: if line has exactly K components, it's a point
        if len(v) == K:
            for k in range(K):
                x[irk,k] = float(v[k])
            try:
                if len(G[str(inm)]) > 0:
                    # only increase rank if point is not isolated
                    irk = irk + 1
            except KeyError:
                print "dg.py:realizationFromFile(): no vtx", str(inm), "possibly isolated"
            inm = inm + 1
            if irk >= n:
                break;
    if irk < n-1:
        print "dg.py:realizationFromFile(): WARNING: rlz in", inputFile, "has <", n, " vtx"
    return x

## read realization from _rlz.dat file 
def realizationFromRlzDatFile(G, inputFile):
    n = G.V.card
    K = NumberOfDimensions
    x = np.zeros((n,K))
    f = open(inputFile, "r")
    offset = 1 # start counting vertices from 1 by default
    idx = 0;
    for line in f:
        v = line.split()
        # assumption: if line has exactly K+1 components, it's a point
        if len(v) == K+1:
            for k in range(1,K+1):
                x[idx,k-1] = float(v[k])
            idx = idx + 1
    return x

############# SECT:BASIC_METHODS ############

## compute the mean scaled distance error of x w.r.t. weights in G
##   if x is a dict, this also works when x is a partial realization of G
##   p is the order of the norm {1,2,inf}
def MDE(G,x,p = 2, scaled=0):
    mde = 0
    for e in G.E:
        tl = e.tail.rank
        hd = e.head.rank
        if type(x) is list:
            xtl = x[tl]
            xhd = x[hd]
        elif type(x) is dict:
            if tl not in x or hd not in x:
                continue
            xtl = x[tl]
            xhd = x[hd]
        else:
            # <type 'numpy.ndarray'>
            xtl = x[tl,:]
            xhd = x[hd,:]
        dij = np.linalg.norm(np.subtract(xtl, xhd), ord=p)
        if e.interval == 0:
            newmde = abs(dij - e.wL)
            if scaled == 1:
                scaling = e.wL
                if scaling < 1:
                    scaling = 1
                newmde = newmde / scaling
        else:
            newmde = max(0, e.wL - dij) + max(0, dij - e.wU)
            if scaled == 1:
                scL = e.wL
                if scL < 1:
                    scL = 1
                scU = e.wU
                if scU < 1:
                    scU = 1
                newmde = max(0, e.wL - dij) / scL + max(0, dij - e.wU) / scU
        mde = mde + newmde
    return mde / len(G.E)

## compute the largest scaled distance error of the w.r.t. weights in G
##   if x is a dict, this also works when x is a partial realization of G
##   p is the order of the norm {1,2,inf}
def LDE(G,x,p = 2, scaled=0):
    lde = 0
    for e in G.E:
        tl = e.tail.rank
        hd = e.head.rank
        if type(x) is list:
            xtl = x[tl]
            xhd = x[hd]
        elif type(x) is dict:
            if tl not in x or hd not in x:
                continue
            xtl = x[tl]
            xhd = x[hd]
        else:
            # <type 'numpy.ndarray'>
            xtl = x[tl,:]
            xhd = x[hd,:]
        dij = np.linalg.norm(np.subtract(xtl, xhd), ord=p)
        if e.interval == 0:
            newlde = abs(dij - e.wL)
            if scaled == 1:
                scaling = e.wL
                if scaling < 1:
                    scaling = 1
                newlde = newlde / scaling
        else:
            newlde = max(0, e.wL - dij) + max(0, dij - e.wU)
            if scaled == 1:
                scL = e.wL
                if scL < 1:
                    scL = 1
                scU = e.wU
                if scU < 1:
                    scU = 1
                newlde = max(0, e.wL - dij) / scL + max(0, dij - e.wU) / scU
        if lde < newlde:
            lde = newlde
    return lde 

## compute the |E|-vector of distance errors of x w.r.t. weights in G
##   if x is a dict, this also works when x is a partial realization of G
##   p is the order of the norm {1,2,inf}
def VDE(G,x, p=2):
    vde = list()
    for e in G.E:
        tl = e.tail.rank
        hd = e.head.rank
        if type(x) is list:
            xtl = x[tl]
            xhd = x[hd]
        elif type(x) is dict:
            if tl not in x or hd not in x:
                continue
            xtl = x[tl]
            xhd = x[hd]
        else:
            # <type 'numpy.ndarray'>
            xtl = x[tl,:]
            xhd = x[hd,:]
        dij = np.linalg.norm(np.subtract(xtl, xhd), ord=p)
        if e.interval == 0:
            vde.append(abs(dij - e.wL))
        else:
            vde.append(max(0, e.wL - dij) + max(0, dij - e.wU))
    return vde 

## verify feasibility of (non-lifted) SDP/DDP (matrix) solution
def MatrixSolutionLDE(G,X, scaled=0):
    n = X.shape[0]
    # check solution
    lde = 0
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        dij = math.sqrt(X[i,i] + X[j,j] - 2*X[i,j] + myZero)
        if e.interval == 0:
            newlde = abs(dij - e.wL)
            if scaled == 1:
                scaling = e.wL
                if scaling < 1:
                    scaling = 1
                newlde = newlde / scaling
        else:
            newlde = max(0,e.wL - dij) + max(0,dij - e.wU)
            if scaled == 1:
                scL = e.wL
                if scL < 1:
                    scL = 1
                scU = e.wU
                if scU < 1:
                    scU = 1
                newlde = max(0,e.wL - dij) / scL + max(0,dij - e.wU) / scU
        if lde < newlde:
            lde = newlde
    return lde

## verify feasibility of (non-lifted) SDP/DDP (matrix) solution
def MatrixSolutionMDE(G,X, scaled=0):
    n = X.shape[0]
    # check solution
    mde = 0
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        dij = math.sqrt(X[i,i] + X[j,j] - 2*X[i,j] + myZero)
        if e.interval == 0:
            newmde = abs(dij - e.wL)
            if scaled == 1:
                scaling = e.wL
                if scaling < 1:
                    scaling = 1
                newmde = newmde / scaling
        else:
            newmde = max(0,e.wL - dij) + max(0,dij - e.wU)
            if scaled == 1:
                scL = e.wL
                if scL < 1:
                    scL = 1
                scU = e.wU
                if scU < 1:
                    scU = 1
                newmde = max(0,e.wL - dij) / scL + max(0,dij - e.wU) / scU
        mde = mde + newmde
    return mde

## return a random realization
def RandomRealization(G, K, M):
    return np.random.uniform(-M,M,(K,G.V.card))

### generate an Achlioptas m x n random projector matrix
def achlioptasProjector(m, n):
    A = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            aij = achlioptas_sample()
            if aij != 0:
                A[i,j] = aij
    return A

### generate a samples from an Achlioptas distribution
def achlioptas_sample():
    ret = 0
    uniformsample = np.random.uniform(0,1)
    if uniformsample < 1.0/6.0:
        ret = -1
    elif uniformsample > 5.0/6.0:
        ret = 1
    else:
        ret = 0
    return ret

## bond angle cosine
def bondAngleCosine(d01, d12, d02):
    costh = (d01**2 + d12**2 - d02**2) / (2*d01*d12)
    # corrections from md-jeep 0.2
    if costh < -1:
        costh = -1
    if costh > 1:
        costh = 1
    return costh

## bond angles cosines given a graph and a DDGP vertex order
##   (this only works for graphs having a strongly DDGP order)
##   theta are the bond angle cos at vertex v spanned by (v-1,v) and (v,v-1),
##     with theta[0] = theta[n-1] = 0 by definition
##   phi are the adjacent angles: arccos phi[v-1] spanned by (v-1,v), (v,v+1)
##     with phi[n-2] = phi[n-1] = 0 by definition
def bondAnglesCosines(G,alpha):
    n = len(alpha)
    theta = [0]
    phi = list()
    hasEdge = True
    for i in range(1,n-1):
        a = G.getEdgeWeight(alpha[i-1],alpha[i])        
        if a == None:
            hasEdge = False
            print "missing edge ({0:d},{1:d})".format(alpha[i-1],alpha[i])
        else:
            a = a[0]
        b = G.getEdgeWeight(alpha[i],alpha[i+1])
        if b == None:
            hasEdge = False
            print "missing edge ({0:d},{1:d})".format(alpha[i],alpha[i+1])
        else:
            b = b[0]
        c = G.getEdgeWeight(alpha[i-1],alpha[i+1])
        if c == None:
            hasEdge = False
            print "missing edge ({0:d},{1:d})".format(alpha[i-1],alpha[i+1])
        else:
            c = c[0]
        if hasEdge:
            thetheta = (a**2 + b**2 - c**2)/(2*a*b)
            if thetheta < -1:
                thetheta = 1
            elif thetheta > 1:
                thetheta = 1
            theta.append(thetheta)
            thephi = (a**2 + c**2 - b**2)/(2*a*c)
            if thephi < -1:
                thephi = 1
            elif thephi > 1:
                thephi = 1
            phi.append(thephi)
    theta.append(0)
    phi.append(0), phi.append(0)
    if len(theta) < n:
        theta = []
        phi = []
    return (theta,phi)

## return the complement (undirected) graph
def complementGraph(G):
    notE = list()
    for uid in range(G.V.card):
        for vid in range(G.V.card):
            u = G.V.byRank[uid]
            v = G.V.byRank[vid]
            if u.rank < v.rank and v not in G[u]:
                notE.append(Edge(u,v, 0, myInf, -1))
    return Graph(E)

## compute ad-bc for a 2x2 matrix [[a b] [c d]]
def det2x2(A):
    return A[0,0]*A[1,1] - A[0,1]*A[1,0]

## dihedral angle cosine (in 3D)
def dihedralAngleCosine3D(d01, d12, d23, d02, d13, d03):
    version = 1

    if version == 0:
        # Leo's version
        costh2 = bondAngleCosine(d01, d12, d02)
        discr2 = 1-costh2**2
        sinth2 = math.sqrt(discr2)
        costh3 = bondAngleCosine(d12, d13, d23) #[not d12,d23,d13] is important!
        discr3 = 1-costh3**2
        sinth3 = math.sqrt(discr3)
        print "sin theta2 =", sinth2, "; sin theta3 =", sinth3
        cosomega = (d01**2 + d13**2 - 2*d01*d13*costh2*costh3 -
                    d03**2) / (2*d01*d13*sinth2*sinth3)
    elif version == 1:
        # Antonio's version (MD-jeep 2.0)
        a = d01*d01 + d13*d13 - d03*d03
        a = a / (2.0*d01*d13)
        b = d13*d13 + d12*d12 - d23*d23
        b = b / (2.0*d13*d12)
        c = d01*d01 + d12*d12 - d02*d02
        c = c / (2.0*d01*d12)
        e = 1.0 - b*b
        f = 1.0 - c*c
        if e < 0.0 or f < 0.0:
            print "dg.py:dihedralAngleCosine3D[MD-jeep]: e =", e, "; f =", f
            exit('dg.py:dihedralAngleCosine3D[MD-jeep]: something wrong') 
        e = math.sqrt(e)
        f = math.sqrt(f)
        cosomega = (a - b*c) / (e*f)
        if cosomega < -1.0:
            cosomega = -1.0
        if cosomega >  1.0:
            cosomega =  1.0
            
    return cosomega



## compute the Gram matrix from the distance matrix
def dist2Gram(EDM):
    n = EDM.shape[0]
    J = np.identity(n) - (1.0/n)*np.ones((n,n))
    B = -0.5 * np.dot(J,np.dot(np.square(EDM), J))
    return B

## compute the distance matrix from the Gram matrix
def Gram2dist(B):
    n = B.shape[0]
    EDM = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                EDM[i,j] = math.sqrt(B[i,i] + B[j,j] - 2*B[i,j])
    return EDM

## return the distance matrix of a realization
##   p is the norm order
def distanceMatrix(x, p=2):
    n = len(x[:,0])
    EDM = np.zeros((n,n))
    for u in range(n-1):
        for v in range(u+1,n):
            EDM[u,v] = np.linalg.norm(np.subtract(x[u,:],x[v,:]), ord=p)
            EDM[v,u] = EDM[u,v]
    return EDM


## matrix without the j-th column
def exceptColumn(A, j):
    return np.delete(A,j,1)


## matrix without the i-th row
def exceptRow(A, i):
    return np.delete(A,i,0)


### generate a Gaussian m x n random projector matrix 
def gaussianProjector(m, n):
    return (1.0/math.sqrt(m)) * np.random.normal([0],[[1]],(m,n))


## return the graph induced by the vertex rank subset S with vtx order R
## S is a set of vertex ranks, so R[S[j]] gives the name of the j-th vertex
## if relabel = False, the subgraph has the same vertex labels as the graph
##   otherwise, it's all relabeled from 0 
def inducedGraphWithOrder(G, S, R, relabelflag=False):
    n = G.V.card
    s = len(S)
    if len(R) == 0:
        R = range(n)
    relabel = {}
    if relabelflag == True:
        if len(R) == n or len(R) == s:
            # DDGP order, just invert R on vertices of S
            i = 0
            for r in R:
                if r in S:
                    relabel[r] = i
                    i = i + 1
        elif len(R) > n:
            # DMDGP order, take first position for repeated vertices
            i = 0
            for r in R:
                if r in S and not (r in relabel):
                    relabel[r] = i
                    i = i + 1
    indE = list()
    for u in S:
        for v in S:
            if u < v:
                if G.isEdge(R[u],R[v]):
                    e = G.getEdge(R[u], R[v])
                    if relabelflag == False:
                        indE.append(e)
                    else:
                        tl = relabel[R[u]]
                        hd = relabel[R[v]]
                        vtl = Vertex(str(tl), tl)
                        vhd = Vertex(str(hd), hd)
                        e2 = Edge(vtl, vhd, e.wL, e.wU, -1)
                        e2.interval = e.interval
                        indE.append(e2)
    GS = Graph(indE)
    return GS

## return the graph induced by the vertex rank subset S
def inducedGraph(G, S, relabel=False):
    return inducedGraphWithOrder(G, S, [], relabel)
    # indE = list()
    # for u in S:
    #     for v in S:
    #         if u < v:
    #             if G.isEdge(u,v):
    #                 indE.append(G.getEdge(u,v))
    # return Graph(indE)


## return True if list of vertex ranks C is a clique in G
def isClique(C,G):
    ret = True
    for u in C:
        for v in C:
            if u < v:
                if not G.isEdge(u,v):
                    ret = False
                    break
        if ret == False:
            break
    return ret

## return True if G is a complete graph
def isComplete(G):
    ret = True
    n = G.V.card
    for u in range(n-1):
        for v in range(u+1,n):
            if not G.isEdge(u,v):
                ret = False
                break
        if ret == False:
            break
    return ret

## is given matrix PSD?
def isPSD(A):
  eigvals, eigvects = np.linalg.eigh(A)
  return np.all(eigvals > -myZero)

## pad an np.array with zeros up to length K
def lift(a, K):
    t = K - len(a)
    if t > 0:
        b = np.pad(a, (0,t), "constant")
    else:
        b = a
    return b

## pad a list with zeros up to length K
def liftList(a, K):
    t = K - len(a)
    if t > 0:
        a.extend([0]*t)
    return a

## get the maximum edge weight from G
def maxEdgeWeight(G):
    return partialDistanceMatrix(G).max()

## count near-zeros of a matrix
def numNearZeros(A):
    ret = 0
    (m,n) = A.shape
    for i in range(m):
        for j in range(n):
            if abs(A[i,j]) <= myZero:
                ret = ret + 1
    return ret

## pair up |det| and numNearZeros
def nzAbsDet(A):
    return [numNearZeros(A), abs(np.linalg.det(A))]

## pair up det2x2 and numNearZeros
def nzDet2x2(A):
    return [numNearZeros(A), det2x2(A)]

## form the partial distance matrix from a graph
##   missing non-diagonal entries have a zero 
def partialDistanceMatrix(G):
    n = G.V.card
    pEDM = np.zeros((n,n))
    for u in range(n):
        for v in G[G.V.byRank[u].name]:
            pEDM[u,G.V.byName[v].rank] = G[G.V.byRank[u].name][v].wL
    return pEDM

## (random) r-diagonal matrices, see [Barvinok DCG 1995, example 4.1]
def rDiagonal(n,r):
    F = np.identity(n)
    for i in range(n):
        for j in range(r):
            if i+j+1 < n:
                rndeps = np.random.uniform(-1.0/r,1.0/r)
                F[i,i+j+1] = rndeps
                F[i+j+1,i] = rndeps
    return F

## range without the i-th element
def rangeBut(r, i):
    R = range(r)
    return R[:i] + R[i+1:]

## return corresponding vertex name list for a list of vertex ranks
def rank2name(G, alpha):
    return [G.V.byRank[v].name for v in range(G.V.card)]

## relative distance (closeness) between two vectors
def relativeCloseness(x,y):
    nx = np.linalg.norm(x)
    ny = np.linalg.norm(y)
    sx = (1/nx) * x
    sy = (1/ny) * y
    angle = abs(np.dot(sx,sy))
    anglefraction = angle / math.pi
    normfraction = min(nx,ny) / max(nx,ny)
    return (anglefraction, normfraction)

## random square positive definite matrix
def rndPositiveDefinite(n):
    F = np.identity(n)
    for i in range(n):
        F[i,i] = np.random.uniform(1.0,5.0)
    P = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            P[i,j] = np.random.uniform(1.0,5.0)
    return np.dot(P.T,np.dot(F,P))

## return the squared distance matrix of a realization
def squareDistanceMatrix(x):
    return np.square(distanceMatrix(x))

## get the sum of all the edge weights (in abs val) of G
def sumEdgeWeight(G):
    return sum(e.wL for e in G.E)

## Householder reflection matrix H w.r.t. vector v normal to the refl hyperpl
##   code: https://people.sc.fsu.edu/~jburkardt/py_src/test_mat/householder.py
def HouseholderReflection(v):
    n = len(v)
    H = np.zeros((n,n))
    for i in range(n):
        H[i,i] = 1.0
    vdot = v.dot(v)
    assert vdot > 0
    for i in range(n):
      for j in range(n):
        H[i,j] = H[i,j] - 2.0*v[i]*v[j]/vdot
    return H

## find normal to a hyperplane
def normal2hyperplane(p):
    # perpendicular v to the plane defined by p
    K = len(p)
    #   as in nextK(), choose a pivot h then form system (p_i - p_h)v = 0
    h = K-1
    A = np.array(map(lambda i:np.subtract(p[i],p[h]), rangeBut(K,h)))
    # A is a K-1 x K matrix, which we expect to be of rank K-1
    # we arbitrarily set v_h = 1 by adding a K-th equation e_h v = 1
    newrow = np.zeros(K)
    newrow[h] = 1
    A = np.vstack([A, newrow])
    # the right hand side is e_K
    b = np.zeros(K)
    b[K-1] = 1
    try:
        v = np.linalg.solve(A,b)
    except LinAlgError:
        print "dg.py:normal2hyperplane(): singular matrix"
        exit('abort')
    return v
        
## find best alignment of two realizations in K dimensions
##   code from http://nghiaho.com/?page_id=671
##   returns a rotation matrix R and a translation t
##   semantics: if y is a rotated/translated version of x, then y = R*x + t
## dimensions of x and y must be n x r where r >= K
##   only the first K components will be considered
def align(x,y,K):
    if type(x) is dict:
        x = dict2array(x)
    if type(y) is dict:
        y = dict2array(y)
    if type(x) is list:
        x = np.array(x)
    if type(y) is list:
        y = np.array(y)
    if len(x) != len(y):
        print "dg.py:align(): |x| =", len(x), "!=", len(y), "= |y|"
        exit('abort')        
    n = x.shape[0]
    r = x.shape[1]
    if r < K:
        print "dg.py:align(): r =", r, "<", K, "= K"
        exit('abort')
    elif r > K:
        print "dg.py:align(): WARNING: r =", r, ">", K, "= K"
        x = x[:,0:K]
        y = y[:,0:K]
    centroid_x = np.mean(x, axis=0)
    centroid_y = np.mean(y, axis=0)
    # centre the points
    xx = x - centroid_x 
    yy = y - centroid_y 
    H = np.transpose(xx).dot(yy)
    U,S,Vt = np.linalg.svd(H)
    R = U.dot(Vt).T
    # special reflection case
    if np.linalg.det(R) < 0:
       # reflection detected
       Vt[K-1,:] *= -1
       R = U.dot(Vt).T
    t = -R.dot(centroid_x.T) + centroid_y.T
    return (R,t)


############# SECT:NEXT_VERTEX ##########

# driver to select the correct next() function
#   p = sequence of k distinct k-dimensional (K-2)-affinely independent points 
#   d = sequence of k distance values, d[i] = dist. to (k-i-1)-th predecessor
def next(p,d,k):
    if k == 1:
        Z = next1D(p,d)
    elif k == 2:
        Z = next2D(p,d)
    elif k == 3:
        Z = next3D(p,d)
        #Z = next3Dtrig(p,d)
        #Z = nextK(p,d,3)
        #Z,evals = nextGram(p,d,3)
    else:
        Z = nextK(p,d,k)
        #Z,evals = nextGram(p,d,k)
    return Z

## return next point from previous point in 1D and its distance
def next1D(p, d):
    return [[p[0,0] - d[0]], [p[0,0] + d[0]]]

## return next point from previous point in 2D and their distances
##    derived formula for general case when b_2 != 0
##    and another formula for special case b_2 = 0
def next2D(p, d):
    Z = list()
    b = np.subtract(p[1],p[0])
    if np.linalg.norm(b) <= myZero:
        print "next2D: p_1 ~= p_0, can't proceed"
        return Z
    if b[1] <= myZero:
        b[1] = 0
    d01sq = b[0]**2 + b[1]**2
    #d01 = math.sqrt(d01sq)
    d02sq = d[0]**2
    d12sq = d[1]**2
    delta012 = (d01sq + d02sq - d12sq) / 2.0
    if b[1] > 0:
        # formula for the general case b_2 != 0
        discr = d01sq * d02sq - delta012**2
        if abs(discr) <= myZero:
            # collinear (one solution only)
            discr = 0
        if discr < -myZero:
            if round(discr) == 0:
                # we're almost collinear, make collinear (one solution only)
                discr = 0
            else:
                print "next2D: metric inequality not satisfied"
                exit('abort')
        z1p = (b[0]*delta012 + b[1]*math.sqrt(discr)) / d01sq
        z2p = math.sqrt(d02sq - z1p**2)
        zp = np.array([z1p, z2p])
        Z.append(np.add(zp, p[0]))
        if discr > 0:
            z1m = (b[0]*delta012 - b[1]*math.sqrt(discr)) / d01sq
            z2m = math.sqrt(d02sq - z1m**2)
            zm = np.array([z1m, z2m])
            Z.append(np.add(zm, p[0]))
    else:
        # if b[1]=0, formula above is invalid, use simpler:
        ## LEO160517: this appears to be wrong on equilateral unit triangle
        z1p = delta012 / b[0]
        z2p = math.sqrt(d02sq - z1p**2)
        zp = np.array([z1p, z2p])
        z1m = z1p
        z2m = -math.sqrt(d02sq - z1p**2)
        zm = np.array([z1m, z2m])
        Z = [zp, zm]
    return Z

## return next point from previous point in 3D and their distances
##   formulae too complicated, go for system, but not in full generality,
##   we use 2x2 inverse matrix formula
def next3D(p, d):
    Z = list()
    # we translate p[0] to the origin first
    p10 = np.subtract(p[1],p[0])
    p20 = np.subtract(p[2],p[0])
    #psq = np.array([np.linalg.norm(p[0])**2, np.linalg.norm(p[1])**2, np.linalg.norm(p[2])**2])
    psq = np.array([0, np.linalg.norm(p10)**2, np.linalg.norm(p20)**2])
    dsq = np.array([d[0]**2, d[1]**2, d[2]**2])
    
    # linear system is A x = b
    #   compute and stabilize A
    A = np.array([ [p10[0], p10[1], p10[2]], [p20[0], p20[1], p20[2]] ])
    # stabilize A: set near-integers to integers
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if abs(A[i,j] - round(A[i,j])) <= myZero:
                A[i,j] = round(A[i,j])
    #  compute and stabilize b
    b = 0.5*np.array([psq[1]-psq[0]-dsq[1]+dsq[0], psq[2]-psq[0]-dsq[2]+dsq[0]])
    for i in range(len(b)):
        if abs(b[i] - round(b[i])) <= myZero:
            b[i] = round(b[i])
    # choose nonbasic: aim at largest number of zeros in nonsingular basics
    nzdet = [nzDet2x2(A[:,[1,2]]), nzDet2x2(A[:,[0,2]]), nzDet2x2(A[:,[0,1]])]
    nbasic = -1
    nz = -1
    det = 0
    for i,nzd in enumerate(nzdet):
        if abs(nzd[1]) > myZero:
            if nzd[0] > nz:
                nz = nzd[0]
                det = nzd[1]
                nbasic = i
    if nbasic == -1:
        # no nonsingular 2x2 submatrices
        exit('dg.py:next3D(): no nonsingular 2x2 submatrix in Ax=b')

    # basic column indices
    basics = range(3)
    basics.remove(nbasic)
    B = A[:,basics]
    
    # try and invert the matrix
    try:
        Binv = np.linalg.inv(B)
    except np.linalg.linalg.LinAlgError as err:
        if 'Singular matrix' in err.message:
            print "dg.py:next3D: weird, nonzero determinant but inversion error"
            return Z
        else:
            raise
    # nonbasic column
    N = A[:,nbasic]
    # alpha = B^-1 N, and remove "close to zero" errors
    alpha = np.dot(Binv, N)
    if abs(alpha[0]) <= myZero:
        alpha[0] = 0
    if abs(alpha[1]) <= myZero:
        alpha[1] = 0    
    # beta = B^-1 b, and remove "close to zero" errors
    beta = np.dot(Binv, b)
    if abs(beta[0]) <= myZero:
        beta[0] = 0
    if abs(beta[1]) <= myZero:
        beta[1] = 0
    if np.linalg.norm(alpha) > myZero:
        # alpha!=0, intersect dictionary with quadratic system
        #   find root of la*x^2 - 2*mu*x + nu = 0
        la = alpha[0]**2 + alpha[1]**2 + 1
        mu = alpha[0]*beta[0] + alpha[1]*beta[1] 
        nu = beta[0]**2 + beta[1]**2 - dsq[0]
        discr = mu**2 - la*nu  # formula (-b\pm sqr{b^2-ac})/a
        if abs(discr) <= myZero:
            # coplanar (one solution only)
            discr = 0
        if discr < -myZero:
            if round(discr) == 0:
                # almost coplanar, make coplanar (one solution only)
                discr = 0
            else:
                # distances are wrong
                print "dg.py:next3D: discr =", discr, "(simplex ineqs not satisfied)"
                exit('abort')
        xp = np.zeros(3)
        xp[nbasic] = (mu + math.sqrt(discr)) / la
        xp[basics[0]] = beta[0] - alpha[0]*xp[nbasic]
        xp[basics[1]] = beta[1] - alpha[1]*xp[nbasic]
        xp = np.add(xp,p[0]) # we'd translated p0 to 0
        Z.append(xp)
        if discr > 0:
            # discriminant != 0, get 2 distinct solutions
            xm = np.zeros(3)
            xm[nbasic] = (mu - math.sqrt(discr)) / la
            xm[basics[0]] = beta[0] - alpha[0]*xm[nbasic]
            xm[basics[1]] = beta[1] - alpha[1]*xm[nbasic]
            xm = np.add(xm,p[0]) # we'd translated p0 to 0
            Z.append(xm)
    else:
        # alpha=0, compute x[basics]=beta first (intro_dg:3.3.2)
        xp = np.zeros(3)
        xp[basics] = beta
        # now x[nbasic] 
        xpnb2 = d[0]**2 - beta[0]**2 - beta[1]**2
        if abs(xpnb2) <= myZero:
            # coplanar, one solution only
            xpnb2 = 0
        if xpnb2 < -myZero:
            if round(xpnb2) == 0:
                # almost coplanar, make coplanar (one solution only)
                xpnb2 = 0
            else:
                # distances are wrong
                print "dg.py:next3D: xpnb2 =", xpnb2, "(simplex ineqs not satisfied)"
                exit('abort')
        xp[nbasic] = math.sqrt(xpnb2) 
        xp = np.add(xp,p[0]) # we'd translated p0 to 0
        Z.append(xp)
        if xpnb2 > myZero:
            # discr != 0, get 2 distinct solutions
            xm = np.zeros(3)
            xm[basics] = beta
            xm[nbasic] = -math.sqrt(xpnb2)
            xm = np.add(xm,p[0]) # we'd translated p0 to 0
            Z.append(xm)
    return Z

# next() in 3D using trigonometric computations
#   improvement: costh could be pre-computed
def next3Dtrig(p, d):
    Z = list()
    # compute rotation matrix U
    v1 = np.subtract(p[2],p[1])
    v1normsq = v1[0]**2 + v1[1]**2 + v1[2]**2
    if v1normsq <= myZero:
        print "dg.py:next3Dtrig: p[1] and p[2] too close"
        return Z
    v1norm = math.sqrt(v1normsq)
    xhat = (1/v1norm)*v1
    v2 = np.subtract(p[1],p[0])
    zhat = np.cross(v1,v2)
    zhat = (1/np.linalg.norm(zhat))*zhat
    yhat = np.cross(xhat,zhat)
    yhat = (1/np.linalg.norm(yhat))*yhat
    U = np.array([xhat, yhat, zhat]).T
    # compute and stabilize bond angle
    costh = bondAngleCosine(v1norm, d[2], d[1])
    if costh > 1:
        costh = 1
    if abs(costh) < myZero:
        costh = 0
    sinth = math.sqrt(1 - costh**2)
    # compute dihedral angle
    v2normsq = v2[0]**2 + v2[1]**2 + v2[2]**2
    if v1normsq <= myZero:
        print "dg.py:next3Dtrig: p[0] and p[1] too close"
        return Z
    v2norm = math.sqrt(v2normsq)
    d02 = np.linalg.norm(np.subtract(p[2],p[0]))
    if d02 <= myZero:
        print "dg.py:next3Dtrig: p[0] and p[2] too close"
        return Z
    cosom = dihedralAngleCosine3D(v2norm, v1norm, d[2], d02, d[1], d[0])
    if cosom > 1:
        cosom = 1
    if abs(cosom) < myZero:
        cosom = 0
    sinomarg = 1 - cosom**2
    if abs(sinomarg) <= myZero:
        # coplanar (one solution only)
        sinomarg = 0
    if sinomarg < -myZero:
        if round(sinomarg) == 0:
            # almost coplanar, make coplanar (one solution only)
            sinomarg = 0
        else:
            print "dg.py:next3Dtrig: arg =", sinomarg, "(can't take sqrt)"
            exit('abort')
    sinom = math.sqrt(sinomarg)
    # compute d-vector and new position
    dd = np.array([-d[2]*costh, d[2]*sinth*cosom, d[2]*sinth*sinom])
    xp = np.add(p[2], np.dot(U,dd))
    Z.append(xp)
    if abs(sinom) > myZero:
        dd = np.array([-d[2]*costh, d[2]*sinth*cosom, -d[2]*sinth*sinom])
        xm = np.add(p[2], np.dot(U,dd))
        Z.append(xm)
    return Z

## return next point from previous points in R^K and their distances - any K
def nextK(p,d,K):
    Z = list()
    # subtract eq(h) from all other eqs ||x-p_i||^2=d_i^2
    h = 0  # best choice of j?
    # system Ax=b
    A = 2 * np.array(map(lambda i:np.subtract(p[i],p[h]), rangeBut(K,h)))
    b = map(lambda i:np.linalg.norm(p[i])**2 - np.linalg.norm(p[h])**2 - (d[i]**2 - d[h]**2), rangeBut(K,h))
    # stabilize A,b: set near-ints to int
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if abs(A[i,j]-round(A[i,j])) <= myZero:
                A[i,j] = round(A[i,j])
    for i in range(len(b)):
        if abs(b[i]-round(b[i])) <= myZero:
            b[i] = round(b[i])
    # choose nonbasic: aim at largest number of zeros/det in nonsing. basics
    nzdet = map(lambda j:nzAbsDet(exceptColumn(A,j)), range(K))
    basics = []
    nbasic = -1
    det = 0
    nz = -1
    for j in range(K):
        if abs(nzdet[j][1]) > myZero:
            if nzdet[j][0] > nz:
                nz = nzdet[j][0]
                det = nzdet[j][1]
                nbasic = j
    if nbasic == -1:
        exit("dg.py:nextK(): no nonsingular 2x2 submatrix of A")
    # basic column indices
    basics = rangeBut(3,nbasic)
    B = A[:,basics]
    # try and invert the matrix
    try:
        Binv = np.linalg.inv(B)
    except np.linalg.linalg.LinAlgError as err:
        if 'Singular matrix' in err.message:
            print "nextK: weird, nonzero determinant but inversion error"
            return Z
        else:
            raise
    # nonbasic column
    N = A[:,nbasic]
    # alpha = B^-1 N, and stabilize it
    alpha = np.dot(Binv, N)
    for j in range(K-1):
        if abs(alpha[j]) <= myZero:
            alpha[j] = 0
    # beta = B^-1 b, and stabilize it
    beta = np.dot(Binv, b)
    for j in range(K-1):
        if abs(beta[j]) <= myZero:
            beta[j] = 0
    if np.linalg.norm(alpha) > myZero:
        # alpha!=0, intersect dictionary with one of the quadratic eqns
        #   select the equation index
        h = 0  # improve choice
        #   find root of la*x^2 - 2*mu*x + nu = 0
        la = np.linalg.norm(alpha)**2 + 1
        mu = np.dot(alpha,beta) - np.dot(alpha,p[h][basics]) + p[h][nbasic]
        nu = np.linalg.norm(beta)**2 + np.linalg.norm(p[h])**2 - 2*np.dot(p[h][basics],beta) - d[h]**2
        discr = mu**2 - la*nu  # formula (-b\pm sqr{b^2-ac})/a
        print la, mu, nu, discr
        if abs(discr) <= myZero:
            # coplanar (one solution only)
            discr = 0
        if discr < -myZero:
            if round(discr) == 0:
                # almost coplanar, make coplanar (one solution only)
                discr = 0
            else:
                # distances are wrong
                print "next: distances do not satisfy simplex ineqs"
                exit('abort')
        xp = np.zeros(K)
        xp[nbasic] = (mu + math.sqrt(discr)) / la
        for j in range(K-1):
            xp[basics[j]] = beta[j] - alpha[j]*xp[nbasic]
        Z.append(xp)
        if discr > 0:
            # discriminant != 0, get 2 distinct solutions
            xm = np.zeros(K)
            xm[nbasic] = (mu - math.sqrt(discr)) / la
            for j in range(K-1):
                xm[basics[j]] = beta[j] - alpha[j]*xm[nbasic]
            Z.append(xm)
    else:
        # alpha=0, compute x[basics] first
        xp = np.zeros(K)
        xp[basics] = beta
        # now x[nbasic] in fn of x[basic] through ||x-p_0||^2=d0^2
        la = 1
        mu = -2*p[0][nbasic]
        nu = np.linalg.norm(beta)**2 - 2*np.dot(p[0][basics],beta) + np.linalg.norm(p[0])**2 - d[0]**2
        discr = mu**2 - la*nu  # formula (-b\pm sqr{b^2-ac})/a
        if abs(discr) <= myZero:
            # coplanar (one solution only
            discr = 0
        if discr < -myZero:
            if round(discr) == 0:
                # almost coplanar, make coplanar (one solution only)
                discr = 0
            else:
                # distances are wrong
                print "next: distances do not satisfy simplex ineqs"
                exit('abort')
        xp[nbasic] = math.sqrt(discr)
        Z.append(xp)
        if discr > myZero:
            # discr > 0, get 2 distinct solutions
            xm = np.zeros(K)
            xm[basics] = beta
            xm[nbasic] = -xp[nbasic]
            Z.append(xm)
    return Z

# this version of next is as follows:
#   slice the partial distance matrix to distances between predecessors in p
#   transform this (complete) distance submatrix into a Gram matrix
#   perform PCA on K largest components, get eigenvectors
#   align with given p's
# differently from other next functions, it returns a pair (Z, evals)
#   use evals[K:] to evaluate realization error
def nextGram(p,d,K):
    Z = list()
    n = len(p)+1
    # square distance matrix (n+1 x n+1)
    D2 = np.zeros((n,n))
    for i in range(n-2):
        for j in range(i+1,n-1):
            dd = np.subtract(p[i],p[j])
            D2[i,j] = dd.dot(dd)
            D2[j,i] = D2[i,j]
    for i in range(n-1):
        D2[i,n-1] = d[i]*d[i]
        D2[n-1,i] = d[i]*d[i]
    # Gram matrix (n+1 x n+1)
    J = np.identity(n) - (1.0/n)*np.ones((n,n))
    B = -0.5 * np.dot(J, np.dot(D2, J))
    # realization over n first columns of the factor of B: x is n by (n-1)
    (evals,evecs) = np.linalg.eigh(B)
    # deal with some numerical errors
    evl = np.zeros(n)
    rank = 0
    for i,ev in enumerate(evals):
        if abs(ev) < myZero:
            evl[i] = 0
        elif ev < -myZero:
            evl[i] = 0
        else:
            evl[i] = ev
            rank += 1
    sqrootdiag = np.eye(n)
    for i in range(n):
        sqrootdiag[i,i] = math.sqrt(evl[i])
    X = np.fliplr(evecs.dot(sqrootdiag))
    evecs = np.fliplr(evecs)
    # if rank == K, there should be no error (only a small floating point one)
    # if rank > K, the metric is not realizable in K dims, we return closest
    #   approximation by ignoring >0 evl beyond K-th and the evl vector
    # if rank < K, the metric is realizable in (rank<K) dimensions and
    #   there's no reflection
    x = X[:,0:K]
    # work out the alignment of x[0], ..., x[n-2] with p[0], ..., p[n-2]
    # since x comes from a factor obtained in decreasing eigenvalue order,
    #   the first rank columns cannot be zero; hence we can project p onto
    #   R^rank and align in that space; and then project the result onto
    #   the subspace spanned by p
    if rank < K:
        pp = np.array(p)
        # first rank components of x are expressed in the standard basis
        #   so projection = restriction to first rank components
        pr = pp[:,0:rank]
        # restore original norm in projected plane
        for i in range(rank):
            prnm = np.linalg.norm(pr[i])
            if prnm > myZero:
                pr[i] = (pr[i] / np.linalg.norm(pr[i])) * np.linalg.norm(pp[i])
        # alignment in rank dimensions
        R,t = align(x[0:n-1,0:rank], pr, rank)
        zp = (R.dot(x[n-1,0:rank].T) + t).T
        # now we should re-project onto the subspace spanned by p, 
        #   ...but I'm exhausted. Pad with zeros and hope for the best
        z = np.hstack((zp, np.zeros(K-rank)))
    else:
        R,t = align(x[0:n-1], p, K)
        # apply alignment to last vector of x
        z = (R.dot(x[n-1,:].T) + t).T
        
    for j in range(K):
        if abs(z[j]) < myZero:
            z[j] = 0
    Z.append(z)
    if rank >= K and len(p) == K:
        ## reflect z across hyperplane def'd by p (if len(p) == K)
        #  1. construct reflection matrix
        v = normal2hyperplane(p)
        H = HouseholderReflection(v)
        #  2. translate z by p_0
        zbar = z - p[0]
        #  3. reflect and translate back
        z = H.dot(zbar) + p[0]
        for j in range(K):
            if abs(z[j]) < myZero:
                z[j] = 0
        Z.append(z)
    return (Z, np.flip(evals, 0))

############# SECT:BRANCH-AND-PRUNE #########

@total_ordering
class KeyDict(object):
    def __init__(self, key, dct):
        self.key = key
        self.dct = dct
    def __lt__(self, other):
        return self.key < other.key
    def __eq__(self, other):
        return self.key == other.key
    def __repr__(self):
        return '{0.__class__.__name__}(key={0.key}, dct={0.dct})'.format(self)

class BPNode:
    def __init__(self, x):
        # partial realization so far
        self.x = x
        # its length
        self.level = len(x)
        # did this node pass through the deleted queue?
        self.delq = False
        # pruning error that got the node in the deleted queue
        self.prunerr = 0

# used to turn the realizations represented by dict to array representation
def dict2array(x):
    return np.array([np.array(x[i]) for i in x])
def array2dict(x):
    return {i:x for i,x in enumerate(list(x))}

# norm2 error between two rlzs (sum of norm2 err on vectors) 
def rlzSumDist(x,y):
    nx = len(x)
    ny = len(y)
    if nx != ny:
        print "np.py:rlzDist(x,y): len(x) =", nx, "!= len(y) =", ny
        exit('abort')
    if type(x) is dict:
        x = dict2array(x)
    if type(y) is dict:
        y = dict2array(y)
    return sum([np.linalg.norm(np.subtract(x[i],y[i])) for i in range(nx)])
        
## the (iterative) Branch-and-Prune
class BPOptions:
    def __init__(self):
        # the pruning distance (PrDst) tolerance
        self.bpeps = 1e-3
        # False: PrDst is absolute error
        self.relerr = True
        # False: bpeps can be exceeded
        self.strictPrune = True
        # 1=some printout, 2=full printout
        self.verbose = 0
        # True: call local NLP solver to "straighten" 
        self.refine = False
        #   (over subproblems of size refwinsize)
        self.refwinsize = 10
        # stop after finding numsol sols (0=all)
        self.numsol = 0
        
def BranchPrune(dgp, opt = None):
    if opt is None:
        opt = BPOptions()    

    ## initialization
    bpeps = opt.bpeps
    K = NumberOfDimensions
    dgp.ordCheck
    G = dgp.graph
    n = G.V.card
    R = dgp.order
    Rrk = [G.V.byName[r].rank for r in R]
    n_ord = len(R)
    # limit length of deleted queue to linear (othw exp)
    maxdelqsize = n_ord
    if opt.verbose > 0:
        print "bp(0): instance has", dgp.type, dgp.ordertype
        print "  order (names):", R
        print "  preds (inst rks):"
        for i in range(K,n_ord):
            print "  ", i, dgp.pred[i]
    instype = dgp.type    # either "ddgp" or "dmdgp"
    # map BP_rank to vertex names from instance
    bp2nm = {i:r for i,r in enumerate(R)}
    nm2bp = {r:[] for r in bp2nm.values()}
    for i in bp2nm:
        nm2bp[bp2nm[i]].append(i)
    # map BP_rank to vertex ranks from instance
    bp2rk = {i:G.V.byName[bp2nm[i]].rank for i in bp2nm}
    rk2bp = {v:[] for v in bp2rk.values()}
    for i in bp2rk:
        rk2bp[bp2rk[i]].append(i)

    ## prepare adjacent predecessor structures
    # adjacent predecessors are partitioned into discretization and pruning    
    #   discretization orders are stored in discr_rk, discr_bp
    #     discr_rk: BP_rank -> list of instance_ranks
    #     discr_bp: BP_rank -> list of BP_ranks
    #   pruning orders are stored in prune_rk, prune_bp
    #     prune_rk: BP_rank -> list of instance_ranks
    #     prune_bp: BP_rank -> list of BP_ranks
    if dgp.ordertype == "order":
        # assumption: len(rk2bp[v])==1 for each v
        # strategy for choosing adj preds: closest K
        discr_rk = {}
        for i in bp2rk:
            dpreds = dgp.pred[bp2rk[i]][-K:]
            discr_rk[i] = [u for j,u in enumerate(dpreds)]
        discr_bp = {i:[rk2bp[u][0] for u in discr_rk[i]] for i in bp2rk}
        
        # pruning predecessors
        prune_rk = {i:[] for i in bp2rk}
        prune_bp = {i:[] for i in bp2rk}
        for i in bp2rk:
            for u in G[G.V.byRank[bp2rk[i]].name]:
                j = rk2bp[G.V.byName[u].rank][0]
                if j < i and j not in discr_bp[i]:
                    prune_rk[i].append(bp2rk[j])
                    prune_bp[i].append(j)
    elif dgp.ordertype == "reorder":
        discr_rk = {i:[] for i in bp2rk}
        discr_bp = {i:[] for i in bp2rk}
        for i in bp2rk:
            if i >= K:
                for j in range(K):
                    discr_rk[i].append(bp2rk[i-j-1])
                    discr_bp[i].append(i-j-1)
        # pruning predecessors
        prune_rk = {i:[] for i in bp2rk}
        prune_bp = {i:[] for i in bp2rk}
        for i in bp2rk:
            for u in G[G.V.byRank[bp2rk[i]].name]:
                try:
                    jl = rk2bp[G.V.byName[u].rank]
                except KeyError:
                    print "dg.py:BranchPrune: warning: vtx name", G.V.byName[u].name, "not in order, ignoring"
                    jl = []
                if len(jl) == 1:
                    # pruning predecessors can't be repetition vertices
                    j = jl[0]
                    if j < i-K:
                        prune_rk[i].append(bp2rk[j])
                        prune_bp[i].append(j)
                        
    #### THE BP ALGORITHM USING A PRIORITY QUEUE
    ## BP initialization
    xinit = dict()
    # set of solutions
    X = list()
    # find position of the first K-clique in order
    if not isClique([bp2rk[i] for i in range(K)], G):
        exit("dg.py:BranchPrune(): first K vtx (in order) not a clique")
    if opt.verbose > 0:
        print "bp(0): initial_clique", [bp2nm[i] for i in range(K)], [bp2rk[i] for i in range(K)]
    xK = RealizeClique(G,K)
    for i in range(K):
        xinit[bp2rk[i]] = xK[i]
    # the priority queue
    bpq = []
    bpn = BPNode(xinit)
    heapq.heappush(bpq, KeyDict(0,bpn))
    # auxiliary queue for non-strict pruning (save nodes there)
    delq = []
    iteration = 0

    ## BP loop
    while len(bpq) > 0:

        ## extract node from the queue
        kdct = heapq.heappop(bpq)
        prty = kdct.key
        node = kdct.dct
        # the partial realization to extend
        xbar = node.x  
        if node.delq and not opt.strictPrune:
            # node comes from deleted queue, change bpeps
            if bpeps < node.prunerr:
                if opt.verbose > 0:
                    print "  bp(", level, "): updating bpeps", bpeps, "->", node.prunerr
                bpeps = node.prunerr

        ## define BP level
        level = len(xbar) # we start ranks from 0
        if opt.verbose > 0:
            print "bp(", level, "): starting iteration", iteration
        iteration = iteration + 1

        ## is xbar a whole realization?
        if level == n_ord:
            ## yes, store it
            if opt.verbose > 0:
                print "bp(", level, "): found new realization"
            X.append(xbar)
            if opt.numsol > 0:
                if len(X) >= opt.numsol:
                    break
            continue

        ## check if this is a repeated vertex
        if dgp.ordertype == "reorder" and len(rk2bp[bp2rk[level]]) > 1:
            ## yes, just go on
            if opt.verbose > 1:
                print "bp(", level, "): repeated vtx inst_rk", bp2rk[level], "with name", bp2nm[level], ": pass"
            continue

        ## BP branching
        vtxrk = bp2rk[level]  # instance rank of current vertex
        vtxnm = bp2nm[level]  # name of current vertex
        if opt.verbose > 0:
            print "  bp(",level,"): find posn of vtx: inst_rk", vtxrk, ", name", vtxnm

        ## discretization predecessors 
        discr_p = [xbar[i] for i in discr_rk[level]]
        # check that discretization predecessors are distinct
        discrOK = True
        for i in discr_rk[level]:
            if discrOK == True:
                for j in discr_rk[level]:
                    if i < j:
                        if np.linalg.norm(np.subtract(xbar[i],xbar[j])) < bpeps:
                            print "  bp(", level, "): WARNING: discr preds", i, j, "too close: pruning this node"
                            discrOK = False
                            break
        if not discrOK:
            continue
        
        discr_edge = [G.getEdge(vtxrk,u) for u in discr_rk[level] if u != vtxrk]
        discr_intvdist = [(e.wL, e.wU, e.interval) for e in discr_edge]
        discr_dist = []
        # verify intervals
        for d in discr_intvdist:
            if not d[2]:
                discr_dist.append(d[0])
            else:
                # discretization distance is an interval
                print "dg.py:BranchPrune(): discretization dist can't be an interval (yet)"
                exit('abort')
        if opt.verbose > 1:
            print "  bp(", level, "): preds", discr_bp[level], discr_rk[level], [bp2nm[j] for j in discr_bp[level]]
            for i in discr_rk[level]:
                print "    dsc_rk", i, ":", list(xbar[i])
            print "  bp(", level, "): discretization distances", discr_dist

        ## find points for current vertex (the "next" call)
        Z = next(discr_p, discr_dist, K)
        if opt.verbose > 0:
            print "  bp(", level, "): found", len(Z), "point(s)"
        if opt.verbose > 1:
            print "  bp(", level, "): pts =", str(Z).replace("array", "").replace("[","").replace("]","")
        # check points are OK w.r.t. distances
        verifydist = [[np.linalg.norm(np.subtract(z, discr_p[i])) for i in range(K)] for z in Z]
        for i,d in enumerate(verifydist):
            # check vectors of distances from z and given
            derr = np.linalg.norm(np.subtract(d,discr_dist))
            if derr > bpeps: ## use BP epsilon for discr errors too
                print "dg.py:BranchPrune(",level,"): next() ret'd pt", i, "w/discrdist err", derr
                ### LEO170508 WARNING: remember to uncomment (if commented)!
                exit('abort')

        ## refinement (call continuous method on fixed window of partial rlz)
        if opt.refine == True:
            S = [bp2rk[i] for i in range(level) if i>=level-opt.refwinsize]
            S = list(set(S))
            GS = inducedGraphWithOrder(G, S, Rrk)
            xstart = np.array([xbar[i] for i in S])
            # method here is local NLP solver, but it can change
            xGS = LocalNLPSolver(GS, xstart, method="lbfgs")
            if opt.verbose > 0:
                rlzerr = rlzSumDist(xGS, xstart)
                print "  bp(", level, "): refinement improved by", rlzerr
            for i,s in enumerate(S):
                xbar[s] = xGS[i]
                
        ## BP pruning
        if len(prune_bp[level]) > 0:

            # pruning predecessors
            prune_p = [xbar[i] for i in prune_rk[level]]
            prune_edge = [G.getEdge(vtxrk,i) for i in prune_rk[level]]
            prune_intvdist = [(e.wL, e.wU, e.interval) for e in prune_edge]
            if opt.verbose > 1:
                print "  bp(", level, "): prunepreds", prune_rk[level], [bp2nm[i] for i in prune_bp[level]]
                for i in prune_rk[level]:
                    print "    prn_rk", i, ":", list(xbar[i])
                print "  bp(", level, "): prunedist", [(rk,prune_intvdist[i][0]) for i,rk in enumerate(prune_rk[level])]

            # compute errors w.r.t. pruning distances
            Zidx = range(len(Z))
            Zdel = [0 for j in Zidx]
            Zerr = [0 for j in Zidx]
            for j in Zidx:
                z = Z[j]
                for i,rk in enumerate(prune_rk[level]):
                    # compare new point z with i-th pruning predecessor 
                    pd = prune_intvdist[i]
                    d = np.linalg.norm(np.subtract(z,xbar[rk]))
                    if opt.verbose > 1:
                        print "  bp( {0:d} ): ||x{1:d} - z{2:d}|| = {3:.3f}".format(level, rk, j, d)

                    # compute distance errors
                    if pd[2]: # interval distance
                        if opt.relerr:
                            derr = max(0,pd[0]-d)/max(pd[0],1)+max(0,d-pd[1])/max(pd[1],1) #relative
                        else:
                            derr = max(0,pd[0]-d)+max(0,d-pd[1]) #absolute
                    else:     # exact distance
                        if opt.relerr:
                            derr = abs(d-pd[0])/max(pd[0],1) #relative
                        else:
                            derr = abs(d-pd[0]) #absolute
                    Zerr[j] += derr

                    # check if distances warrant pruning
                    if derr > bpeps:
                        if opt.verbose > 0:
                            print "  bp(", level, "): strict pruning on sol", j,
                            if opt.relerr:
                                print "by rel. err",
                            else:
                                print "by abs. err",
                            print derr, ">", bpeps                            
                        # YES: record j on the deleted list
                        Zdel[j] = 1  
                        if opt.strictPrune:
                            # strict pruning: we can stop with one wrong dist
                            break
                        else:
                            # if pruning is not strict, need all distances
                            pass

            # see if everything was deleted
            if not opt.strictPrune and sum(Zdel) == len(Z):
                # slack pruning and everything deleted
                for j in Zidx:
                    if len(delq) < maxdelqsize:
                        Zerr[j] = Zerr[j]/level # average PrDst error
                        x = {i:xbar[i] for i in xbar} # extend xbar by Z[j]
                        x[bp2rk[level]] = Z[j]
                        if opt.verbose > 1:
                            print "  bp(", level, "): adding next level node to delq, prty", 1/Zerr[j]
                        # add node to deleted queue
                        bpd = BPNode(x)
                        bpd.prunerr = Zerr[j]
                        bpd.delq = True
                        newprty = 1/Zerr[j] # get smallest error first
                        heapq.heappush(delq, KeyDict(newprty, bpd))

            # delete pruned points from Z
            Z = [z for i,z in enumerate(Z) if Zdel[i] == 0]

        else:
            # no pruning distances
            if opt.verbose > 0:
                print "  bp(", level, "): no pruning distances"
            
        # make new nodes in the priority queue
        serprty = 1    # this privileges "depth first" (serprty is increase)
        topprty = 1e5  # this is the maximum priority if LDE is "zero"
        for z in Z:
            # extend x=(xbar,z)
            x = {i:xbar[i] for i in xbar} 
            x[bp2rk[level]] = z
            # add node to queue
            bpn = BPNode(x)
            # priority is 1/LDE
            newprty = LDE(G,x)
            if newprty < myZero:
                newprty = topprty
            else:
                newprty = 1 / newprty 
            newprty += serprty
            heapq.heappush(bpq, KeyDict(newprty, bpn))
            serprty += 1
        if opt.verbose > 0:
            print "  bp(", level, "): adding", len(Z), "node(s) to queue"

        # deal with no solutions, empty bpq and slack pruning
        if not opt.strictPrune and len(bpq) == 0 and len(X) == 0:
            # add deleted nodes in the queue
            delkdct = heapq.heappop(delq)
            heapq.heappush(bpq, delkdct)            
    return (X, bp2nm, bp2rk)


## this version of the Branch-and-Prune is ONLY for nextGram()
##   nextGram() can take all predecessors, so there's no explicit pruning
##   nor for partitioning predecessors into discretization and pruning
def BranchPruneGram(dgp, opt = None):
    if opt is None:
        opt = BPOptions()
    ## initialization
    bpeps = opt.bpeps
    K = NumberOfDimensions
    dgp.ordCheck
    G = dgp.graph
    n = G.V.card
    R = dgp.order
    Rrk = [G.V.byName[r].rank for r in R]
    n_ord = len(R)
    # limit length of deleted queue to linear (othw exponential)
    maxdelqsize = n_ord
    if opt.verbose > 0:
        print "bpG(0): instance has", dgp.type, dgp.ordertype
        print "  order (names):", R
        print "  preds (inst rks):"
        for i in range(K,n_ord):
            print "  ", i, dgp.pred[i]
    instype = dgp.type    # either "ddgp" or "dmdgp"
    # map BP_rank to vertex names from instance
    bp2nm = {i:r for i,r in enumerate(R)}
    nm2bp = {r:[] for r in bp2nm.values()}
    for i in bp2nm:
        nm2bp[bp2nm[i]].append(i)
    # map BP_rank to vertex ranks from instance
    bp2rk = {i:G.V.byName[bp2nm[i]].rank for i in bp2nm}
    rk2bp = {v:[] for v in bp2rk.values()}
    for i in bp2rk:
        rk2bp[bp2rk[i]].append(i)

    ## prepare adjacent predecessor structures
    #    pred_rk: BP_rank -> list of instance_ranks
    #    pred_bp: BP_rank -> list of BP_ranks
    if dgp.ordertype == "order":
        # assumption: len(rk2bp[v])==1 for each v
        pred_rk = {}
        for i in bp2rk:
            dpreds = dgp.pred[bp2rk[i]]
            pred_rk[i] = [u for j,u in enumerate(dpreds)]
        pred_bp = {i:[rk2bp[u][0] for u in pred_rk[i]] for i in bp2rk}        
    elif dgp.ordertype == "reorder":
        # repetition orders, take contiguous K first, then "pruning"
        pred_rk = {i:[] for i in bp2rk}
        pred_bp = {i:[] for i in bp2rk}
        for i in bp2rk:
            if i >= K:
                for j in range(K):
                    pred_rk[i].append(bp2rk[i-j-1])
                    pred_bp[i].append(i-j-1)
        for i in bp2rk:
            for u in G[G.V.byRank[bp2rk[i]].name]:
                try:
                    jl = rk2bp[G.V.byName[u].rank]
                except KeyError:
                    print "dg.py:BranchPrune: warning: vtx name", G.V.byName[u].name, "not in order, ignoring"
                    jl = []
                if len(jl) == 1:
                    # pruning predecessors can't be repetition vertices
                    j = jl[0]
                    if j < i-K:
                        pred_rk[i].append(bp2rk[j])
                        pred_bp[i].append(j)
                        
    #### THE BP ALGORITHM USING A PRIORITY QUEUE AND nextGram()
    ## BP initialization
    xinit = dict()
    # set of solutions
    X = list()
    # find position of the first K-clique in order
    if not isClique([bp2rk[i] for i in range(K)], G):
        exit("dg.py:BranchPrune(): first K vtx (in order) not a clique")
    if opt.verbose > 0:
        print "bpG(0): initial_clique", [bp2nm[i] for i in range(K)], [bp2rk[i] for i in range(K)]
    xK = RealizeClique(G,K)
    for i in range(K):
        xinit[bp2rk[i]] = xK[i]
    # the priority queue
    bpq = []
    bpn = BPNode(xinit)
    heapq.heappush(bpq, KeyDict(0,bpn))
    # auxiliary queue for non-strict pruning (save nodes there)
    delq = []
    iteration = 0

    ## BP loop
    while len(bpq) > 0:

        ## extract node from the queue
        kdct = heapq.heappop(bpq)
        prty = kdct.key
        node = kdct.dct
        xbar = node.x
        if node.delq and not opt.strictPrune:
            # node comes from deleted queue, update epsilon
            if bpeps < node.prunerr:
                if opt.verbose > 0:
                    print "  bpG(", level, "): updating bpeps", bpeps, "->", node.prunerr
                bpeps = node.prunerr

        ## define BP level
        level = len(xbar) # we start ranks from 0
        if opt.verbose > 0:
            print "bpG(", level, "): starting iteration", iteration
        iteration = iteration + 1

        ## is xbar a whole realization?
        if level == n_ord:
            ## yes, store it
            if opt.verbose > 0:
                print "bpG(", level, "): found new realization"
            X.append(xbar)
            if opt.numsol > 0:
                if len(X) >= opt.numsol:
                    break
            continue

        ## check if this is a repeated vertex
        if dgp.ordertype == "reorder" and len(rk2bp[bp2rk[level]]) > 1:
            ## yes, just go on
            if opt.verbose > 1:
                print "bpG(", level, "): repeated vtx inst_rk", bp2rk[level], "with name", bp2nm[level], ": pass"
            continue

        ## BP branching
        vtxrk = bp2rk[level]  # instance rank of current vertex
        vtxnm = bp2nm[level]  # name of current vertex
        if opt.verbose > 0:
            print "  bpG(",level,"): find posn of vtx: inst_rk", vtxrk, ", name", vtxnm

        ## adjacent predecessors 
        pred_p = [xbar[i] for i in pred_rk[level]]
        # the "if u != vtxrk" is for re-orders: count preds only once
        pred_edge = [G.getEdge(vtxrk,u) for u in pred_rk[level] if u != vtxrk]
        pred_intvdist = [(e.wL, e.wU, e.interval) for e in pred_edge]
        pred_dist = []
        # verify intervals
        for d in pred_intvdist:
            if not d[2]:
                pred_dist.append(d[0])
            else:
                # predecessor distance is an interval
                print "dg.py:BranchPrune(): distances can't be intervals (yet)"
                exit('abort')
        if opt.verbose > 1:
            print "  bpG(", level, "): preds", pred_bp[level], pred_rk[level], [bp2nm[j] for j in pred_bp[level]]
            for i in pred_rk[level]:
                print "    pred_rk", i, ":", list(xbar[i])
            print "  bpG(", level, "): predecessor distances", pred_dist

        ## find points for current vertex (ONLY nextGram() allowed here)
        (Z,evl) = nextGram(pred_p, pred_dist, K)
        if opt.verbose > 0:
            print "  bpG(", level, "): found", len(Z), "point(s)"
        if opt.verbose > 1:
            print "  bpG(", level, "): pts =", str(Z).replace("array", "").replace("[","").replace("]","")
            print "            evals:", evl
            
        # "pruning" (simply check computed distances fit given ones)
        Zidx = range(len(Z))
        Zdel = [0 for j in Zidx]
        Zerr = [0 for j in Zidx]
        for j in Zidx:
            z = Z[j]
            for i,rk in enumerate(pred_rk[level]):
                # compare new point z with i-th adjacent predecessor 
                pd = pred_intvdist[i]
                d = np.linalg.norm(np.subtract(z,xbar[rk]))
                if opt.verbose > 1:
                    print "  bpG( {0:d} ): ||x{1:d} - z{2:d}|| = {3:.3f}".format(level, rk, j, d)

                # compute distance errors
                if pd[2]: # interval distance
                    if opt.relerr:
                        derr = max(0,pd[0]-d)/max(pd[0],1)+max(0,d-pd[1])/max(pd[1],1) #relative
                    else:
                        derr = max(0,pd[0]-d)+max(0,d-pd[1]) #absolute
                else:     # exact distance
                    if opt.relerr:
                        derr = abs(d-pd[0])/max(pd[0],1) #relative
                    else:
                        derr = abs(d-pd[0]) #absolute
                Zerr[j] += derr

                # check if distance errors warrant pruning
                if derr > bpeps:
                    if opt.verbose > 0:
                        print "  bpG(", level, "): strict pruning on sol", j,
                        if opt.relerr:
                            print "by rel. err",
                        else:
                            print "by abs. err",
                        print derr, ">", bpeps                            
                    # YES: record j on the deleted list
                    Zdel[j] = 1  
                    if opt.strictPrune:
                        # strict pruning: we can stop with one wrong dist
                        break
                    else:
                        # if pruning is not strict, need all distances
                        pass

        # see if "pruning" deleted everything
        if not opt.strictPrune and sum(Zdel) == len(Z):
            # slack pruning and everything deleted
            for j in Zidx:
                if len(delq) < maxdelqsize:
                    Zerr[j] = Zerr[j]/level # average PrDst error
                    x = {i:xbar[i] for i in xbar} # extend xbar by Z[j]
                    x[bp2rk[level]] = Z[j]
                    if opt.verbose > 1:
                        print "  bpG(", level, "): adding next level node to delq, prty", 1/Zerr[j]
                    # add node to deleted queue
                    bpd = BPNode(x)
                    bpd.prunerr = Zerr[j]
                    bpd.delq = True
                    newprty = 1/Zerr[j] # get smallest error first
                    heapq.heappush(delq, KeyDict(newprty, bpd))

        # delete pruned points from Z
        Z = [z for i,z in enumerate(Z) if Zdel[i] == 0]
                
        ## refinement (call continuous method on fixed window of partial rlz)
        if opt.refine == True:
            S = [bp2rk[i] for i in range(level) if i >= level-opt.refwinsize]
            S = list(set(S))
            GS = inducedGraphWithOrder(G, S, Rrk)
            xstart = np.array([xbar[i] for i in S])
            # method here is local NLP solver, but it can change
            xGS = LocalNLPSolver(GS, xstart, method="lbfgs")
            if opt.verbose > 0:
                rlzerr = rlzSumDist(xGS, xstart)
                print "  bp(", level, "): refinement improved by", rlzerr
            for i,s in enumerate(S):
                xbar[s] = xGS[i]
                
        # make new nodes in the priority queue
        serprty = 1    # this privileges "depth first" (serprty is increase)
        topprty = 1e5  # this is the maximum priority if LDE is "zero"
        for z in Z:
            # extend x=(xbar,z)
            x = {i:xbar[i] for i in xbar} 
            x[bp2rk[level]] = z
            # add node to queue
            bpn = BPNode(x)
            # priority is 1/LDE
            newprty = LDE(G,x)
            if newprty < myZero:
                newprty = topprty
            else:
                newprty = 1 / newprty 
            newprty += serprty
            heapq.heappush(bpq, KeyDict(newprty, bpn))
            serprty += 1
        if opt.verbose > 0:
            print "  bpG(", level, "): adding", len(Z), "node(s) to queue"

        # deal with no solutions, empty bpq and slack pruning
        if not opt.strictPrune and len(bpq) == 0 and len(X) == 0:
            # add deleted nodes in the queue
            delkdct = heapq.heappop(delq)
            heapq.heappush(bpq, delkdct)            
    return (X, bp2nm, bp2rk)


############# SECT:COMBINATORIAL #########

## find a DDGP order in G
##   this returns the tuple (nm_ord,rk_ord,YES/NO,map,inv_map)
##   where "nm_ord" is the DDGP order as a list of vertex names
##         "rk_ord" is the DDGP order as a list of vertex ranks
##         "YES/NO" returns the True/False as an anser to the decision problem
##         "map" maps original vtx rk order to DDGP order equivalent
##         "inv_map" is the inverse map
##   if boolean return flag is NO then nameTOPorder is empty but rankTOPorder
##      gets as far as the alg could extract DDGP order off last clique
def DDGPOrder(G, K = NumberOfDimensions):
    n = G.V.card
    for C in itertools.combinations(range(n), K):
        alpha = list()
        lC = list(C)
        sC = set(C)
        if isClique(C,G):
            ordFromThisClique = True
            a = dict()
            for v in range(n):
                a[v] = 0
            # start the order with the clique
            alpha = alpha + lC
            W = set(range(n)) - sC
            for v in W:
                Nv = set([G.V.byName[u].rank for u in G[G.V.byRank[v].name]])
                a[v] = len(Nv & sC)
            while len(W) > 0:
                (v,adjPr) = max([(u,a[u]) for u in W],
                                key=operator.itemgetter(1))
                if adjPr < K:
                    ordFromThisClique = False
                    break
                alpha.append(v)
                for u in G[G.V.byRank[v].name]:
                    urk = G.V.byName[u].rank
                    if urk in W:
                        a[urk] = a[urk] + 1
                W.remove(v)
            if ordFromThisClique:
                nat2top = dict()
                for i in range(n):
                    nat2top[i] = alpha[i]
                top2nat = dict()
                for i in range(n):
                    top2nat[alpha[i]] = i
                R = [G.V.byRank[v].name for v in alpha]
                return (R, alpha, ordFromThisClique, nat2top, top2nat)
    return ([], alpha, ordFromThisClique, None, None)


## Floyd-Warshall:
##   return EDM with missing edge weights replaced by shortest path lengths 
def FloydWarshall(G):
    n = G.V.card
    M = maxEdgeWeight(G)
    EDM = (M+1)*np.ones((n,n))
    # initialize with existing weights
    for u in range(n):
        # zero diagonal
        EDM[u,u] = 0
        for v in G[G.V.byRank[u].name]:
            EDM[u,G.V.byName[v].rank] = G[G.V.byRank[u].name][v].wL
    for u in range(n):
        for v in range(n):
            for w in range(n):
                d = EDM[u,w] + EDM[w,v]
                if EDM[u,v] > d:
                    EDM[u,v] = d
    return EDM

## Floyd-Warshall:
##   return EDM with missing edge weights replaced by shortest path lengths
##   this version has a default big M = (sum edge weights)
def FloydWarshall1(G):
    n = G.V.card
    M = sumEdgeWeight(G)
    EDM = (M+1)*np.ones((n,n))
    # initialize with existing weights
    for u in range(n):
        # zero diagonal
        EDM[u,u] = 0
        for v in G[G.V.byRank[u].name]:
            EDM[u,G.V.byName[v].rank] = G[G.V.byRank[u].name][v].wL
    for u in range(n):
        for v in range(n):
            for w in range(n):
                d = EDM[u,w] + EDM[w,v]
                if EDM[u,v] > d:
                    EDM[u,v] = d
    return EDM

## Floyd-Warshall:
##   return EDM with missing edge weights replaced by shortest path lengths
##   this version uses the fortran-compiled code in dgopt.so
def FloydWarshall2(G):
    n = G.V.card
    M = maxEdgeWeight(G)
    pEDM = (M+1)*np.ones((n,n))
    # initialize with existing weights
    for u in range(n):
        # zero diagonal
        pEDM[u,u] = 0
        for v in G[G.V.byRank[u].name]:
            pEDM[u, G.V.byName[v].rank] = G[G.V.byRank[u].name][v].wL
    flatpEDM = pEDM.flatten()
    nn = n*n
    dgopt.floydwarshall2(n, nn, flatpEDM)
    EDM = np.reshape(flatpEDM, (n,n), order='C')
    return EDM

## the IsoMap method
def IsoMap(G):
    print "dg.py:IsoMap(): Floyd-Warshall (compiled version)..."
    A = FloydWarshall2(G)
    print "dg.py:IsoMap(): EDM to Gram..."
    B = dist2Gram(A)
    print "dg.py:IsoMap(): Principal Component Analysis..."
    x = PCA(B, NumberOfDimensions)
    return x

## the IsoMap method with a full distance matrix
def IsoMapEDM(A):
    B = dist2Gram(A)
    x = PCA(B, NumberOfDimensions)
    return x

## this version of IsoMap skips PCA and returns (x,rank)
def IsoMapFullRank(G):
    A = FloydWarshall(G)
    B = dist2Gram(A)
    x = PCA(B)
    return (x,rk)

## the IsoMap heuristic mixes classic IsoMap, JLL and local NLP
##   as of 160224, this is a failure: bottleneck is NLP local solution
##   also, JLL doesn't seem to be doing much good
def IsoMapHeuristic(G, locnlp = "ipopt"):
    jlleps = 0.1
    jllconst = 1
    n = G.V.card
    jllk = int(round(jllconst * math.log(n) / (jlleps*jlleps)))
    if n <= jllk:
        print "IsoMapHeuristic: |V| <", jllk, " too small, calling IsoMap(G)"
        return IsoMap(G)
    K = NumberOfDimensions
    # MDS and keep all nonzero eigenvalues >= jlleps * max_eigenvalue
    A = FloydWarshall(G)
    B = dist2Gram(A)
    (evals,evecs) = np.linalg.eigh(B)
    eig_threshold = jlleps * math.sqrt(max(evals)) # fraction of sqr(max eig)
    evals[evals < eig_threshold] = 0
    rk = sum(1 for eig in evals if eig > myZero)
    X = MDS(B,eig_threshold)
    if rk < K:
        # if not enough columns, pad with zeros
        X = np.append(X, np.zeros((n,K-rk)), axis=1)
    # project cols of x using JLL
    print "IsoMapHeuristic: projecting from", n, " to", jllk
    # sample a jllk x n random projector matrix
    T = gaussianProjector(jllk,rk)
    TX = np.dot(T,X)
    # now adjust TX to the original distances (round of NLP)
    #TX1 = LocalNLPSolver(G, TX, {}, method=locnlp, maxitn = 10) ## need speed
    #TX1 = StochasticProximityEmbedding(G, TX, p=0.01)
    TX1 = TX ## no optimization at all
    print LDE(G,TX1)
    # get the Gram matrix of TX1
    TB = np.dot(TX1.T,TX1)
    # perform PCA
    x1 = PCA(TB,K)
    #x = LocalNLPSolver(G, x1, {}, method=locnlp, maxitn = 5) ## need speed
    #x = StochasticProximityEmbedding(G, x1, p=0.01)
    x = x1
    return x

## a very simply multistart metaheuristic based on IPOpt
def MultiStart(G, iterations, method="ipopt"):
    # size
    K = NumberOfDimensions
    n = G.V.card
    # bounds [-M,M]
    A = FloydWarshall(G)
    M = max(A.flatten())
    # obj fun
    ldestar = sum(A.flatten())    
    # MS
    print "MultiStart(", method, "):"
    print "itn ldeinit lde ldestar"
    for t in range(iterations):
        xinit = RandomRealization(G,K,M)
        ldeinit = LDE(G,xinit)
        x = LocalNLPSolver(G,xinit,{}, method)
        lde = LDE(G,x)
        if lde < ldestar:
            ldestar = lde
            xstar = x
        print t, ldeinit, lde, ldestar
    return xstar

## factor a square matrix (new version 170308)
def factor(A):
    n = A.shape[0]
    (evals,evecs) = np.linalg.eigh(A)
    evals[evals < 0] = 0  # closest SDP matrix
    X = evecs #np.transpose(evecs)
    sqrootdiag = np.eye(n)
    for i in range(n):
        sqrootdiag[i,i] = math.sqrt(evals[i])
    X = X.dot(sqrootdiag)
    return np.fliplr(X)

## perform classic Multidimensional scaling
def MDS(B, eps = myZero):
    n = B.shape[0]
    x = factor(B)
    (evals,evecs) = np.linalg.eigh(x)
    K = len(evals[evals > eps])
    if K < n:
        # only first K columns
        x = x[:,0:K]
    return x

## perform PCA (new version 170309)
def PCA(B, K = "None"):
    x = factor(B)
    n = B.shape[0]
    if type(K) == types.StringType:
        K = n
    if K < n:
        # only first K columns
        x = x[:,0:K]
    return x

# Frechet's embedding in given dimension:
#   complete distance matrix to find Frechet in n dims, then
#   use MILP to find the best choice of K columns off n
def Frechet(G, K, timeLim = 120):
    t0 = time.time()
    n = G.V.card
    A = FloydWarshall(G) # complete the graph
    P = pic.Problem()
    y = P.add_variable('y', n, vtype = 'binary')
    s = P.add_variable('s', (n,n))
    w = []
    for k in range(n):
        w.append(P.add_variable('w[%d]'%k, (n,n), vtype = 'binary'))
    # fix non-edge variables and add bounds
    for i in range(n):
        for j in range(n):
            if G.isEdge(i,j):
                P.add_constraint(s[i,j] > 0.0)
                P.add_constraint(s[i,j] < np.linalg.norm(A[i]-A[j], np.inf))
            else:
                P.add_constraint(s[i,j] == 0.0)
    P.set_objective('min', sum(s[i,j] for i in range(n) for j in range(n) if G.isEdge(i,j)))
    # constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        P.add_constraint(sum(w[k][i,j] for k in range(n)) > 1.0)
        for k in range(n):
            P.add_constraint(abs(A[i,k]-A[j,k])*y[k] - s[i,j] < A[i,j])
            P.add_constraint(abs(A[i,k]-A[j,k])*y[k] + s[i,j] > A[i,j]*w[k][i,j])
    P.add_constraint(sum(y[k] for k in range(n)) == K)
    t1 = time.time()
    #P.write_to_file('dgp1.lp')
    sol = P.solve(solver = 'cplex', timelimit = timeLim, cplex_params = {'mip.display':2})
    yval = [i for (i,v) in enumerate(list(y.value)) if v == 1]
    #print y.value, yval
    x = A[:,yval]
    t2 = time.time()
    cpu = t2 - t0
    print "Frechet: LDEinf=", LDE(G,x,np.inf), ", MDEinf=", MDE(G,x,np.inf), ", CPU=", cpu
    return x

# Matousek's embeddings (randomized Frechet) for 1-norm, followed by PCA
#   this heuristic produces crap results (Leo160706)
def Matousek(G, K):
    t0 = time.time()
    A = FloydWarshall(G) # complete the graph
    Amax = A.max()
    n = G.V.card
    D = math.log(n)
    m = int(round(math.pow(n,(2.0/D)*math.log(n))))
    p = math.pow(n,(-2.0/D))
    S = dict()  # the sets S_ij in Racke's lecture notes
    d = dict()  # distances from each v to sets S_ij
    jmax = int(math.ceil(D/2.0))
    dimMax = jmax*m
    y = list()
    for j in range(jmax):
        for i in range(m):
            S[(i,j)] = set()
            for h in range(n):
                d[(h,i,j)] = Amax
                if np.random.uniform() < math.pow(p,j):
                    S[(i,j)].add(h)
    # define the embedding using the S sets
    for v in range(n):
        yv = np.zeros(dimMax)
        for j in range(jmax):
            for i in range(m):
                dim = j*jmax + i
                for h in S[(i,j)]:
                    if A[v,h] < d[(v,i,j)] - myZero:
                        d[(v,i,j)] = A[v,h]
                yv[dim] = d[(v,i,j)]
        y.append(yv)
    y = np.array(y)
    Y = dist2Gram(distanceMatrix(y, p=1))
    y = PCA(Y,K)
    cpu = time.time() - t0
    print "Matousek: LDE1 =", LDE(G,y,1), ", MDE1 =", MDE(G,y,1), " CPU =", cpu
    return y

# Frechet embedding, then JLL projection
def RandomProjection(G):
    A = FloydWarshall(G)  # complete the graph
    # trivially, for all i <= n, A[i] is the Frechet embedding of vector i
    # now project using JLL
    jlleps = 0.15
    jllconst = 1.8
    n = G.V.card
    jllk = int(round(jllconst * math.log(n) / (jlleps*jlleps)))
    #jllk = NumberOfDimensions
    K = NumberOfDimensions
    if n <= jllk:
        print "RandomProjection: |V| <", jllk, " too small"
        return np.zeros((n,K))
    else:
        print "RandomProjection: projecting from", n, "to", jllk
    T = gaussianProjector(jllk, n)
    TA = np.dot(T, A)
    x = PCA(np.dot(TA.T, TA), K)
    return x

## realize a K-clique in R^K (arbitrarily choose second vertex for each pair)
#  (assumes that the first K vertices in G are a clique)
def RealizeClique(G,K):
    x = []
    # position 1st vertex
    x = np.zeros((K,K))
    for k in range(1,K):
        eStr = [G.adj[G.V.byRank[h].name][G.V.byRank[k].name] for h in range(k)]
        d = [e.wL for e in eStr]
        Z = next(x[0:k], d, k)
        x[k] = lift(Z[1], K) # pad with zeros requires x initialized at [0]
    return x


## stochastic proximity embedding algorithm
##   x is the starting realization (dim x numpoints matrix)
##   number of cycles at each iteration (i.e. cardinality of
##   of subset of terms in stochastic gradient descent)
##   is p|E|
def StochasticProximityEmbedding(G, x, p = 0.01):
    learning = 1.0
    m = len(G.E)
    C = int(round(p * m))
    learnStep = p
    spemultfact = 16
    while learning > 0:
        for c in range(C):
            ij = int(round(np.random.uniform(0,m)))
            if ij >= m:
                ij = m-1
            e = G.E[ij]
            i = e.tail.rank
            j = e.head.rank
            nij = spemultfact * (np.dot(x[i,:],x[i,:]) + np.dot(x[i,:],x[i,:]) - 2*np.dot(x[i,:],x[j,:]))
            dij = e.wL
            x[i,:] = x[i,:] + (learning*(dij - nij) / nij)*np.subtract(x[i,:],x[j,:])
            x[j,:] = x[j,:] + (learning*(dij - nij) / nij)*np.subtract(x[j,:],x[i,:])
        learning = learning - learnStep
    return x

## verify whether the vertex order of the given instance is a K-DMDGP order 
def isDMDGPInstance(dgp,K):
    G = dgp.graph
    R = dgp.order
    n = G.V.card
    ret = False
    S = range(K)  # set of ranks
    # verify initial clique
    if not isComplete(inducedGraphWithOrder(G,S,R)):
        return ret
    # verify each vertex has at K adjacent (contiguous) predecessors
    for i,v in enumerate(R):
        if i >= K:
            ret = True
            for u in R[i-K:i-1]:
                if not (G.isEdgeByName(u,v) or u == v):
                    print "dg.py:isDMDGPInstance():", u, v, "not edge or =="
                    return False
    return True

## verify whether the vertex order of the graph G is a K-DDGP order 
def isDDGPOrder(G,K):
    n = G.V.card
    ret = False
    S = range(K)
    # verify initial clique
    if not isComplete(inducedGraph(G, S)):
        return ret
    # verify each vertex has at least K adjacent predecessors
    for v in range(K,n):
        Nv = set([G.V.byName[u].rank for u in G[G.V.byRank[v].name]])
        if len(Nv & set(range(v))) < K:
            return ret
    return True

## spanning tree based constructive heuristic
##   grow a breadth-first spanning tree based on largest degree first
##   realize vertices in the BFS order of the tree, each vertex randomly
##   on the sphere S(parent, edge length)
def SpanningTreeEmbedding(G,K):
    n = G.V.card
    K = NumberOfDimensions
    x = np.zeros((n,K))
    marked = np.zeros(n)  # marked vertices are already in the tree
    degSeq = [(len(G[G.V.byRank[u].name]),u) for u in range(n)]
    maxDeg = max(degSeq)
    deg = dict()
    for (d,v) in degSeq:
        deg[v] = d
    q = Q.PriorityQueue()
    q.put((-maxDeg[0],maxDeg[1]))
    R = set()
    R.add(maxDeg[1])
    while not q.empty():
        v = q.get()[1]
        for wname in G[G.V.byRank[v].name]:
            w = G.V.byName[wname].rank
            if w not in R:
                if G.isEdge(v,w):
                    dist = G.getEdgeWeight(v,w)
                elif G.isEdge(w,v):
                    dist = G.getEdgeWeight(w,v)
                rnddir = np.random.uniform(-1,1,K)
                rnddir /= np.linalg.norm(rnddir)
                rnddir *= 0.5*(dist[0]+dist[1])
                x[w,:] = np.add(x[v,:],rnddir)
                R.add(w)
                q.put((-(deg[w]),w))
    return x

############# SECT:FORMULATION_BASED #############

## solve a DGP by SDP/DDLP, run Barvinok's naive algorithm, then locnlpsolve
def Barvinok(G, method = "sdp", ddp_itn = 5, locnlp="lbfgs", sdpFormulation=1):
    n = G.V.card
    K = NumberOfDimensions
    if method == "sdp":
        sol = SDPModelAndSolve(G, sdpFormulation)
    elif method == "ddp":
        sol = DDLPFullRank(G, sdpFormulation)[0]
    elif method == "iterddp":
        sol = IterativeDDLP(G, ddp_itn, sdpFormulation)
    else:
        print "Barvinok: method should be in {\"sdp\", \"ddp\"}"
        exit
    # clean up solutions from slightly negative eigenvalues and factor it
    t1 = time.time()
    X = factor(sol)
    t2 = time.time()
    print "Barvinok: time to factor =", t2-t1
    # Barvinok's naive algorithm:
    barvlde = myInf
    xs = np.zeros((n,K))
    for cnt in range(BarvinokIterations):
        # sample from nK-multivariate normal random vector
        y = (1/math.sqrt(K))*np.random.multivariate_normal(np.zeros(n*K), np.identity(n*K))
        #y = np.random.uniform(-1,1,nK)
        y = np.reshape(y, (K, n))
        # x = y X
        x1 = np.transpose(np.dot(y,X)) # not X.T, factor() output already tr
        x1lde = LDE(G, x1)
        if x1lde < barvlde:
            xs = x1
            barvlde = x1lde
    t1 = time.time()
    x = LocalNLPSolver(G, xs, {}, locnlp)
    t2 = time.time()
    print "Barvinok: time to call local solver =", t2-t1
    print "Barvinok: after Barvinok's alg, LDE =", barvlde
    print "Barvinok: after locnlp", locnlp, ", LDE =", LDE(G,x)    
    return x

## Barvinok's naive algorithm v.2.0 in [B., Course in Convexity, p.242]
def BarvinokBook(G, method = "sdp", ddp_itn = 5, locnlp="lbfgs", sdpFormulation=1):
    n = G.V.card
    K = NumberOfDimensions
    if method == "sdp":
        sol = SDPModelAndSolve(G, sdpFormulation)
    elif method == "ddp":
        sol = DDLPFullRank(G, sdpFormulation)[0]
    elif method == "iterddp":
        sol = IterativeDDLP(G, ddp_itn, sdpFormulation)
    else:
        print "BarvinokBook: method should be in {\"sdp\", \"ddp\"}"
        exit
    # clean up solutions from slightly negative eigenvalues and factor it
    t1 = time.time()
    X = factor(sol)
    t2 = time.time()
    print "BarvinokBook: time to factor =", t2-t1
    # Barvinok's naive algorithm v.2 (p.242 of his book):
    barvlde = myInf
    xs = np.zeros((n,K))
    Y0 = np.zeros(n)
    ym = np.random.multivariate_normal(np.zeros(n*1), np.identity(n*1))
    ym = np.reshape(ym, (1, n))
    for cnt in range(BarvinokIterations):
        # sample from n1-multivariate normal random vector
        y = np.random.multivariate_normal(np.zeros(n*1), np.identity(n*1))
        y = np.reshape(y, (1, n))
        Y0 = Y0 + y.T.dot(ym)
    Y0 = Y0 / BarvinokIterations
    X0 = X.T.dot(Y0.dot(X))
    xs = PCA(X0, NumberOfDimensions)
    t1 = time.time()
    barvlde = LDE(G, xs)
    x = LocalNLPSolver(G, xs, {}, locnlp)
    t2 = time.time()
    print "BarvinokBook: time to call local solver =", t2-t1
    print "BarvinokBook: after Barvinok's alg v.2, LDE =", barvlde
    print "BarvinokBook: after locnlp", locnlp, ", LDE =", LDE(G,x)    
    return x

## Barvinok's low-rank SDP solution [Barvinok DCG 1995]
def BarvinokSDP(G,r):
    n = G.V.card
    #return SDPEDMCP(G,rDiagonal(G.V.card,r))
    return SDPwithObjective(G,rDiagonal(G.V.card,r))

## solve a DGP by diag dominant approximaxation and
##   return the K-dimensional solution projection slice
def DDLP(G, objtype = 1):
    n = G.V.card
    K = NumberOfDimensions
    U = np.identity(n+K)
    sol = DDLPModelAndSolve(G,U, objtype)
    x = PCA(sol, K)
    return x

## this version of DDLP skips the PCA and returns (x,rank)
def DDLPFullRank(G, objtype = 1):
    n = G.V.card
    K = NumberOfDimensions
    U = np.identity(n+K)
    sol = DDLPModelAndSolve(G, U, objtype)
    x = PCA(sol)
    rk = np.linalg.matrix_rank(x)
    return (x, rk)

## build and solve the DDP innder approx of the SDP rlx of DGP in PICOS
##   objtype == 1: Leo
##   objtype == 2: Ye
##   objtype == 3: Barvinok '95
## LEO170615: there appears to be a PICOS bug --
##   iterddp yields increasingly small obj fun values when maximizing
##   matlab code ddprealize/iterddp behaves differently (larger obj fun vals)
def DDLPModelAndSolve(G, U = "None", objtype = 1):
    n = G.V.card
    K = NumberOfDimensions
    t0 = time.time()
    P = pic.Problem()
    if type(U) == types.StringType:
        U = np.identity(K+n)
    elif type(U) is np.ndarray:
        if U.shape[0] != U.shape[1]:
            print "dg.py:DDLPModelAndSolve(): U not square, shape is", U.shape
            exit('abort')
        if U.shape[0] not in [n, n+K]:
            print "dg.py:DDLPModelAndSolve(): U is not n x n or n+K x n+K"
            exit('abort')
        if U.shape[0] == n:
            IK = np.identity(K)
            Zn = np.zeros((K,n))
            U = np.bmat([[IK, Zn], [U, Zn.T]])
    else:
        print "dg.py:DDLPModelAndSolve(): type(U) =", type(U), "!= np.ndarray"
        exit('abort')
        
    IdK = cvx.matrix(np.identity(K))

    # decision variables
    X = P.add_variable('X', (n,n), vtype='symmetric')
    x = P.add_variable('x', (n,K))
    Z = P.add_variable('Z', (n+K,n+K), vtype='symmetric')
    T = P.add_variable('T', (n+K,n+K), vtype='symmetric')
    #S = P.add_variable('S', (n,n), vtype='symmetric')

    ## objective function
    if objtype == 1:
        ## Leo (push-and-pull)
        P.set_objective('max',pic.sum([X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] for e in G.E], 'e', 'G.E'))
        ## min slacks
        #P.set_objective('min', sum(T[i,j] for i in range(n+K) for j in range(n+K)) + sum(S[i,j] for i in range(n) for j in range(n)))
    elif objtype == 3:
        ## Barvinok 1995
        P.set_objective('min', ( cvx.matrix(rDiagonal(n,K)) | X ))
        #P.set_objective('min', ( cvx.matrix(rndPositiveDefinite(n)) | X ))
        
    ## constraints

    ## zero slacks (not necessary)
    #for i in range(n):
    #    for j in range(n):
    #        if not G.isEdge(i,j):
    #            P.add_constraint(S[i,j] == 0)
    
    # distance constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        if objtype == 1:
            ## use with "push-and-pull" objective
            if e.interval == 0:
                P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] <= e.wL**2)
            else:
                P.add_constraint(e.wL**2<=X[i,i]+X[j,j]-2*X[i,j]<=e.wU**2)
        else:
            ## sandwiching
            #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 >= -S[i,j])
            #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 <= S[i,j])
            ## only allow a positive slack -- enough to counter DDP infeas
            #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 == S[i,j])
            if e.interval == 1:
                P.add_constraint(e.wL**2 <= X[i,i] + X[j,j] - 2*X[i,j] <= e.wU**2)
            else:
                P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] == e.wL**2)

    # Z is diagonally dominant
    for i in range(0,n+K):
        P.add_constraint(sum(T[i,j] for j in range(n+K) if j != i) <= Z[i,i])

    # Schur complement = U Z U for some DD Z
    Up = cvx.matrix(np.matrix(U))
    UpT = cvx.matrix(np.matrix(U.T))
    P.add_constraint( (( IdK & x.T ) // ( x & X )) == UpT * Z * Up )
    P.add_constraint( Z <=  T )
    P.add_constraint( Z >= -T )
    #P.add_constraint( T >= 0 )
    #P.add_constraint( S >= 0 )
    t1 = time.time()
    # solve and get solution
    try:
        P.solve(solver = 'cplex', cplex_params = {'lpmethod': 4, 'barrier.crossover': -1})
        #P.solve(solver = 'mosek')
        Xval = np.array(X.value)
    except cplex.exceptions.errors.CplexSolverError:
        Xval = rDiagonal(n, K)
    #Yval = np.array(UpT * Z.value * Up)
    print "DDP: PSD(X) =", isPSD(Xval), ", LDE =", MatrixSolutionLDE(G, Xval)
    t2 = time.time()
    print "DDP: modelling time =", t1-t0, "; solving time =", t2-t1
    rksol = np.linalg.matrix_rank(Xval)    
    ddpobj = sum(Xval[e.head.rank,e.head.rank] + Xval[e.tail.rank,e.tail.rank] - 2*Xval[e.head.rank,e.tail.rank] for e in G.E)
    print "DDP: sol. rank = {0:d}, obj fun val = {1:f}".format(rksol, ddpobj);
    return Xval

## iteratively solve DDLP approximations to the classic SDP relaxation
##    this version calls DDLPModelAndSolve() repeatedly
##    note that U here is n x n (not (K+n)x(K+n))
def IterativeDDLP(G, iterations = 5, objtype = 1):
    n = G.V.card
    K = NumberOfDimensions
    U = np.identity(n)
    t0 = time.time()
    for i in range(iterations):
        ## solve the DDP
        Xv = DDLPModelAndSolve(G, U, objtype)
        objval = sum(Xv[e.head.rank,e.head.rank] + Xv[e.tail.rank,e.tail.rank] - 2*Xv[e.head.rank,e.tail.rank] for e in G.E)
        lde = LDE(G,PCA(Xv,K))
        print "IterativeDDLP(itn", i+1, "/", iterations, "): obj =", objval, "rank sol = ", np.linalg.matrix_rank(Xv), "LDE(X) = {0:0.2f}".format(MatrixSolutionLDE(G,Xv)), "LDE(x) = {0:0.2f}".format(lde)
        if lde < myZero:
            print "  realization error in K dim is <", myZero, ", stopping"
            break
        ## prepare next iteration
        ## Cholesky factor
        #Xv += 1e-5*np.eye(n)
        #U = np.linalg.cholesky(Xv)
        ## factor
        U = factor(Xv)
        ## ignore the eigenvalues
        #(evals,U) = np.linalg.eigh(Xv)
    t1 = time.time()
    print "IterativeDDLP: solving time =", t1-t0
    return Xv

## iteratively solve DDLP approximations to the classic SDP relaxation
##    this version is independent of DDLPModelAndSolve(), and attempts
##    to change the existing formulation in memory; it works very
##    poorly because PICOS takes ages to remove_constraint()
def IterativeDDLPWarmStart(G, iterations = 5, objtype = 1):
    n = G.V.card
    K = NumberOfDimensions
    U = np.identity(n+K)
    P = pic.Problem()
    IdK = cvx.matrix(np.identity(K))

    t0 = time.time()
    # decision variables
    X = P.add_variable('X', (n,n), vtype='symmetric')
    x = P.add_variable('x', (n,K))
    Z = P.add_variable('Z', (n+K,n+K), vtype='symmetric')
    T = P.add_variable('T', (n+K,n+K), vtype='symmetric')
    #S = P.add_variable('S', (n,n), vtype='symmetric')

    ## objective function
    ## push-and-pull
    P.set_objective('min',pic.sum([X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] for e in G.E], 'e', 'G.E'))
    ## min slacks S,T
    #P.set_objective('min', sum(T[i,j] for i in range(n+K) for j in range(n+K)) + sum(S[i,j] for i in range(n) for j in range(n)))
    ## min slacks T
    #P.set_objective('min', sum(T[i,j] for i in range(n+K) for j in range(n+K)))
        
    ## constraints

    ## zero slacks (not necessary)
    #for i in range(n):
    #    for j in range(n):
    #        if not G.isEdge(i,j):
    #            P.add_constraint(S[i,j] == 0)
    
    # distance constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        # to use with "push-and-pull" objective
        P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] >= e.wL**2)
        ## sandwiching
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 >= -S[i,j])
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 <= S[i,j])
        ## only allow a positive slack -- enough to counter DDP infeas
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 == S[i,j])

    # Z is diagonally dominant
    for i in range(n+K):
        P.add_constraint(sum(T[i,j] for j in range(n+K) if j != i) <= Z[i,i])
    P.add_constraint( Z <=  T )
    P.add_constraint( Z >= -T )
    P.add_constraint( T >= 0 )
    #P.add_constraint( S >= 0 )

    # Schur complement = U Z U for some DD Z --- put it last
    Up = cvx.matrix(np.matrix(U))
    UpT = cvx.matrix(np.matrix(U.T))
    P.add_constraint( (( IdK & x.T ) // ( x & X )) == UpT * Z * Up )
    
    t1 = time.time()    
    print "IterativeDDLP: modelling time =", t1-t0

    for i in range(iterations):
        ## solve the DDP
        P.solve(solver = 'cplex', cplex_params = {'lpmethod': 4, 'barrier.crossover': -1})
        ## get solution values
        Xv = np.array(X.value)
        Yv = np.array(UpT * Z.value * Up)

        ## compute error
        #lde = MatrixSolutionLDE(G,X)
        lde = LDE(G, PCA(Xv,K))
        ofval = P.obj_value()
        print "IterativeDDLP(itn", i+1, "/", iterations, "): rank sol = ", np.linalg.matrix_rank(Xv), "LDE = {0:0.2f}".format(lde), ", objfun = {0:0.2f}".format(ofval)
        if lde < myZero:
            break

        ## prepare next iteration
        U = factor(Yv)
        P.remove_constraint(len(P.constraints)-1) # Schurcompl constraint last
        Up = cvx.matrix(np.matrix(U))
        UpT = cvx.matrix(np.matrix(U.T))
        P.add_constraint( (( IdK & x.T ) // ( x & X )) == UpT * Z * Up )
        
    print "IterativeDDLP: PSD(X) =", isPSD(Xv), ", LDE =", MatrixSolutionLDE(G, Xv)
    t2 = time.time()
    print "IterativeDDLP: solving time =", t2-t1
    return Xv


## (this version of LocalNLPSolver is from dg_isomap.py --
##  the one from formulations.py is in old/LocalNLPSolver.py)
##
## call to local solver for min sum_uv (||xu-xv||^2 - d_uv^2)^2
##   G is the input graph, xinit is the "starting point", and
##   xfixed is a map vertex->points to be fixed before optimization
##     this allows optimizing only some points in the realization
##
def LocalNLPSolver(G, xinit, xfixed = {}, method="ipopt", maxitn = nlpMaxIter):
    # this class serves the purpose of making member attributes available
    #    to subfunctions (just a tech detail to overcome a Python 2.x limit)
    class local:
        gph = G
        fxd = xfixed # this is a dictionary of lists not a list of lists
                     # initialize to [] if no fixed vertices
        K = NumberOfDimensions #xinit.shape[1]

    # evaluate objective function sum_uv (||x_u-x_v||^2-d_uv^2)^2
    #    warning: x is a flat array, needs to be mapped to a structured
    #       array with indices in V and K
    def eval_fun(x):
        for v in local.fxd:
            for k in range(local.K):
                x[v*local.K + k] = local.fxd[v][k]
        ret = 0
        for e in local.gph.E:
            dij2 = sum((x[e.head.rank*local.K+l] - x[e.tail.rank*local.K+l])**2 for l in range(local.K))
            dij = np.sqrt(dij2)
            if e.interval == 0:
                ret = ret + (dij2 - e.wL**2)**2
                #ret = ret + (dij - 1/e.wL)**2  # inverse -- by Katya
            else:
                ret = ret + (max(0,e.wL**2-dij2))**2 + (max(0,dij2-e.wU**2))**2
                #ret = ret + (max(0,1/e.wL-dij))**2 + (max(0,dij-1/e.wU))**2
        return ret

    # return |E|-vector ||x_u-x_v||^2-d_uv^2
    #    warning: x is a flat array, needs to be mapped to a structured
    #       array with indices in V and K
    def eval_sys(x):
        for v in local.fxd:
            for k in range(local.K):
                x[v*local.K + k] = local.fxd[v][k]
        ret = list()
        for e in local.gph.E:
            dij2 = sum((x[e.head.rank*local.K+l] - x[e.tail.rank*local.K+l])**2 for l in range(local.K))
            dij = np.sqrt(dij2)
            if e.interval == 0:
                ret.append(dij - e.wL)
                #ret.append(dij2 - e.wL**2)
                #ret.append(dij - 1/e.wL)   #inverse -- suggested by Katya
            else:
                ret.append(max(0,e.wL - dij) + max(0,dij - e.wU))
                #ret.append(max(0,e.wL**2 - dij2) + max(0,dij2 - e.wU**2))
                #ret.append(max(0,1/e.wL - dij) + max(0,dij - 1/e.wU))
        return np.array(ret)
    
    # actually do the work
    xflatinit = xinit.flatten()
    if maxitn > 0:
        print "LocalNLPSolver: calling", method, "for", maxitn, "iterations"
    else:
        print "LocalNLPSolver: calling", method

    if method == "ipopt":
        pyipopt.set_loglevel(0)
        # I edited pyipopt/pyipoptpackage/ipoptunconstrained.py
        #   so it accepts nlpMaxIter as a last argument,
        #   see ipoptunconstrained.py in current dir
        if maxitn > 0:
            sol = pyipopt.fmin_unconstrained(eval_fun, xflatinit, functools.partial(eval_grad, eval_fun), functools.partial(eval_hess, eval_fun), maxiter=maxitn)
        else:
            sol = pyipopt.fmin_unconstrained(eval_fun, xflatinit, functools.partial(eval_grad, eval_fun), functools.partial(eval_hess, eval_fun))
        X = sol[0]
    elif method == "lbfgs":
        if maxitn > 0:
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), method='L-BFGS-B', options={'disp':True, 'maxiter': maxitn})
        else:
#            sol = optimize.fmin_l_bfgs_b(eval_fun, xflatinit, fprime=functools.partial(eval_grad, eval_fun)); X = sol[0]
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), method='L-BFGS-B')
        X = sol.x
    elif method == "cg":
        if maxitn > 0:
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), method='CG', options={'disp':True, 'maxiter':maxitn})
        else:
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), method='CG')
        X = sol.x
    elif method == "newtoncg":
        if maxitn > 0:
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), hess=functools.partial(eval_hess, eval_fun), method='Newton-CG', options={'disp':True, 'maxiter':maxitn})
        else:
            sol = optimize.minimize(eval_fun, xflatinit, jac=functools.partial(eval_grad, eval_fun), hess=functools.partial(eval_hess, eval_fun), method='Newton-CG')
        X = sol.x
    elif method == "powell":
        if maxitn > 0:
            sol = optimize.minimize(eval_fun, xflatinit, method='Powell', options={'disp':True, 'maxiter':maxitn})
        else:
            sol = optimize.minimize(eval_fun, xflatinit, method='Powell')
        X = sol.x
    elif method == "neldermead":
        if maxitn > 0:
            sol = optimize.minimize(eval_fun, xflatinit, method='Nelder-Mead', options={'disp': True, 'maxiter':maxitn})
        else:
            sol = optimize.minimize(eval_fun, xflatinit, method='Nelder-Mead')
        X = sol.x
    elif method == "leastsq":
        if maxitn > 0:
            sol = optimize.leastsq(eval_sys, xflatinit, ftol = 1e-2, xtol = 1e-2, maxfev = maxitn)
        else:
            sol = optimize.leastsq(eval_sys, xflatinit, ftol = 1e-2, xtol = 1e-2)
        X = sol[0]
    elif method == "dgsol" or method == "dgopt":
        # this is different from the rest - calling the optimizer in dgsol
        KK = local.K
        m = len(G.E)
        n = G.V.card
        nK = KK*n
        fOpt = 0.0
        if maxitn > 0:
            (X,fopt) = dgopt.dgopt(KK, n, m, nK, xflatinit, 
                                   [G.E[i].tail.rank for i in range(m)],
                                   [G.E[i].head.rank for i in range(m)],
                                   [G.E[i].wL for i in range(m)],
                                   [G.E[i].wU for i in range(m)],
                                   [0 for i in range(max(2*n,m))],
                                   [0.0 for i in range(22*nK + 21)], maxitn)
        else:
            (X,fopt) = dgopt.dgopt(KK, n, m, nK, xflatinit, 
                                   [G.E[i].tail.rank for i in range(m)],
                                   [G.E[i].head.rank for i in range(m)],
                                   [G.E[i].wL for i in range(m)],
                                   [G.E[i].wU for i in range(m)],
                                   [0 for i in range(max(2*n,m))],
                                   [0.0 for i in range(22*nK + 21)], 1e+3)
    # elif method == "nlopt":
    #     nK = G.V.card*local.K
    #     #opt = nlopt.opt(nlopt.LD_SLSQP, nK)
    #     #opt = nlopt.opt(nlopt.LD_MMA, nK)
    #     #opt = nlopt.opt(nlopt.LD_LBFGS, nK)
    #     opt = nlopt.opt(nlopt.LD_TNEWTON_PRECOND_RESTART, nK)
    #     def objfun(x, grad):
    #         grad = np.zeros(nK)
    #         f = eval_fun(x)
    #         x = algopy.UTPM.init_jacobian(x)
    #         grad[:] = algopy.UTPM.extract_jacobian(eval_fun(x))
    #         return f
    #     opt.set_min_objective(objfun)
    #     varbound = sum([e.wU for e in G.E])
    #     opt.set_lower_bounds(-varbound)
    #     opt.set_upper_bounds(varbound)
    #     X = opt.optimize(xflatinit)
    #     fopt = opt.last_optimum_value()
    #     result = opt.last_optimize_result()

    elif method == "root":
        ## not accepted by 'lm' method
        # def jacfun(x):
        #     K = local.K
        #     E = local.gph.E
        #     m = len(E)
        #     grad = np.zeros(m*K)
        #     eidx = 0
        #     for e in E:
        #         i = e.tail.rank
        #         j = e.head.rank
        #         xi = x[K*i:K*(i+1)-1]
        #         xj = x[K*j:K*(j+1)-1]
        #         grad[K*eidx:K*(eidx+1)-1] = np.subtract(xi,xj)
        #         eidx = eidx + 1
        #     return grad
        if maxitn > 0:
            sol = optimize.root(eval_sys, xflatinit, method='lm', options={'ftol':1e-2, 'xtol':1e-3, 'maxiter':maxitn})
        else:
            sol = optimize.root(eval_sys, xflatinit, method='lm', options={'ftol':1e-2, 'xtol':1e-3})
        X = sol.x
    else:
        print "LocalNLPSolver: method", method, " not known"
        sys.exit(115)
    xstar = np.reshape(X, (local.gph.V.card, local.K), order='C')
    return xstar
# used to compute the gradient
def eval_grad(f, theta):
    theta = algopy.UTPM.init_jacobian(theta)
    return algopy.UTPM.extract_jacobian(f(theta))
# used to compute the Hessian
def eval_hess(f, theta):
    theta = algopy.UTPM.init_hessian(theta)
    return algopy.UTPM.extract_hessian(len(theta), f(theta))    


## solve a DGP by SDP and return the K-dimensional solution projection slice
def SDP(G, sdpFormulation = 1):
    sol = SDPModelAndSolve(G, sdpFormulation)
    x = PCA(sol, NumberOfDimensions)
    return x


## EDMCP SDP feasibility relaxation with given objective function min F.Y
def SDPEDMCP(G,F):
    n = G.V.card
    K = NumberOfDimensions
    P = pic.Problem()
    IdK = cvx.matrix(np.identity(K))
    # decision variables
    Y = P.add_variable('Y', (n,n), vtype='symmetric')
    # objective function
    F1 = cvx.matrix(F)
    P.set_objective('min', F1|Y)
    # constraints
    for e in G.E:
        if e.interval == 1:
            P.add_constraint(e.wL**2 <= Y[e.head.rank,e.head.rank] +
                             Y[e.tail.rank,e.tail.rank] -
                             2*Y[e.head.rank,e.tail.rank] <= e.wU**2)
        else:
            P.add_constraint(Y[e.head.rank,e.head.rank] +
                             Y[e.tail.rank,e.tail.rank] -
                             2*Y[e.head.rank,e.tail.rank] == e.wL**2)
    P.add_constraint(Y >> 0)
    # solve and get solution
    if (solverTolerance >= 0):
        P.solve(tol = solverTolerance, mosek_params={'intpnt_tol_dfeas' : solverTolerance})
    else:
        P.solve()
    sol = np.array(Y.value)
    rk = np.linalg.matrix_rank(sol)
    return (sol, rk)

## EDMCP DDP feasibility relaxation with given objective function min F.Y
#### TO DO (working here)
def DDPEDMCP(G,F):
    n = G.V.card
    K = NumberOfDimensions
    P = pic.Problem()
    IdK = cvx.matrix(np.identity(K))
    # decision variables
    X = P.add_variable('Y', (n,n), vtype='symmetric')
    T = P.add_variable('T', (n,n), vtype='symmetric')
    S = P.add_variable('S', (n,n), vtype='symmetric')
    # objective function
    F1 = cvx.matrix(F)
    P.set_objective('min', (F1)|X)
    # distance constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] == e.wL**2)
        ## sandwiching
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 >= -S[i,j])
        #P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 <= S[i,j])
        ## only allow a positive slack -- enough to counter DDP infeas
        if e.interval == 1:
            P.add_constraint(e.wL**2 - S[i,j] <= X[i,i] + X[j,j] - 2*X[i,j] <= e.wU**2 + S[i,j])
        else:
            P.add_constraint(X[i,i] + X[j,j] - 2*X[i,j] - e.wL**2 == S[i,j])
    # Z is diagonally dominant
    for i in range(n):
        P.add_constraint(sum(T[i,j] for j in range(n) if j != i) <= X[i,i])
    P.add_constraint(X <=  T)
    P.add_constraint(X >= -T)
    P.add_constraint(T >= 0)
    P.add_constraint(S >= 0)
    #P.add_constraint(S <= 100)
    # solve and get solution
    t0 = time.time()
    if (solverTolerance >= 0):
        P.solve(tol = solverTolerance, mosek_params={'intpnt_tol_dfeas' : solverTolerance})
    else:
        P.solve()
    cpu = time.time() - t0
    sol = np.array(X.value)
    rk = np.linalg.matrix_rank(sol)
    return (sol, rk, cpu)

## this version of SDP skips the PCA and returns (x,rank)
def SDPFullRank(G, sdpFormulation = 1):
    sol = SDPModelAndSolve(G, sdpFormulation)
    x = PCA(sol)
    rk = np.linalg.matrix_rank(x)
    return (x, rk)

## build and solve the classic SDP relaxation of the DGP in PICOS
##   objtype == 1: Leo
##   objtype == 2: Ye
##   objtype == 3: Barvinok
##   objtype == 4: arbitrary matrix F, stored in objmatrix
##   objtype == 5: minimize largest eigenvalue of X
##   objtype == 6: feasibility only (zero matrix in objective)
## returns the PSD matrix out of the SDP solver
def SDPModelAndSolve(G, objtype = 1, objmatrix = None):
    n = G.V.card
    K = NumberOfDimensions
    t0 = time.time()
    P = pic.Problem()
    IdK = cvx.matrix(np.identity(K))
    # decision variables
    X = P.add_variable('X', (n,n), vtype='symmetric')
    x = P.add_variable('x', (n,K))
    # objective function
    if objtype == 1:
        ## Leo "push-and-pull" 
        P.set_objective('max',pic.sum([X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] for e in G.E], 'e', 'G.E'))
        #P.set_objective('min',pic.sum([X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] for e in G.E], 'e', 'G.E'))
    elif objtype == 2:
        ## Ye
        P.set_objective('min', pic.tracepow(X,1))
    elif objtype == 3:
        ## Barvinok 1995
        P.set_objective('min', ( cvx.matrix(rDiagonal(n,K)) | X ))
        #P.set_objective('min', ( cvx.matrix(rndPositiveDefinite(n)) | X ))
    elif objtype == 4:
        ## arbitrary multiplier matrix in obj fun
        P.set_objective('min', ( cvx.matrix(objmatrix) | X ))
    elif objtype == 5:
        ## largest eigenvalue of X (used with Barvinok's heuristic)
        tmin = P.add_variable('t', (1,1))
        P.set_objective('min', tmin)
        P.add_constraint(pic.sum_k_largest_lambda(X,1) <= tmin)
    #else:
        ## zero matrix in obj fun (feasibility only)
    # constraints
    for e in G.E:
        if objtype == 1:
            ## Leo (push-and-pull): ineq constraints
            if e.interval == 0:
                P.add_constraint(X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] <= e.wL*e.wL)
            else:
                P.add_constraint(e.wL**2 <= X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] <= e.wU**2)
            
        else:
            ## any other: == constraints
            if e.interval == 1:
                P.add_constraint(e.wL**2 <= X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] <= e.wU**2)
            else:
                P.add_constraint(X[e.head.rank,e.head.rank] + X[e.tail.rank,e.tail.rank] - 2*X[e.head.rank,e.tail.rank] == e.wL**2)
    P.add_constraint(( ( IdK & x.T ) // ( x & X ) ) >> 0)
    #P.add_constraint(X >> 0)
    t1 = time.time()
    # solve and get solution
    P.solve()
    t2 = time.time()
    print "SDPModelAndSolve: modelling time =", t1-t0, "; solving time =", t2-t1    
    sol = np.array(X.value)
    rksol = np.linalg.matrix_rank(factor(sol))
    print "SDP(sdpFormulation={0:d}): solution rank = {1:d}".format(objtype, rksol);
    return sol

## SDP feasibility formulation with given objective function min F.Y
def SDPwithObjective(G,F):
    n = G.V.card
    K = NumberOfDimensions
    P = pic.Problem()
    IdK = cvx.matrix(np.identity(K))
    # decision variables
    X = P.add_variable('X', (n,n), vtype='symmetric')
    x = P.add_variable('x', (K,n))
    # objective function
    F1 = cvx.matrix(F)
    P.set_objective('min', F1|X)
    # constraints
    for e in G.E:
        if e.interval == 1:
            P.add_constraint(e.wL**2 <= X[e.head.rank,e.head.rank] +
                             X[e.tail.rank,e.tail.rank] -
                             2*X[e.head.rank,e.tail.rank] <= e.wU**2)
        else:
            P.add_constraint(X[e.head.rank,e.head.rank] +
                             X[e.tail.rank,e.tail.rank] -
                             2*X[e.head.rank,e.tail.rank] == e.wL**2)
    #P.add_constraint(( ( IdK & x ) // ( x.T & X ) ) >> 0)
    P.add_constraint(X >> 0)
    # solve and get solution
    if solverTolerance >= 0:
        P.solve(tol = solverTolerance, mosek_params={'intpnt_tol_dfeas' : solverTolerance})
    else:
        P.solve()
    sol = np.array(X.value)
    rk = np.linalg.matrix_rank(sol)
    return (sol, rk)


## 160224 this is not working out as expected on large-scale instances
##
## within IsoMapHeuristic:
##    TX1 = StochasticNLPSolver(G, TX, {})
##
## stochastic descent for min sum_uv (||xu-xv||^2 - d_uv^2)^2
##   G is the input graph, xinit is the "starting point", and
##   xfixed is a map vertex->points to be fixed before optimization
##     this allows optimizing only some points in the realization
def StochasticNLPSolver(G, xinit, xfixed = {}, p = 0.01, close = 0.001):
    # this class serves the purpose of making member attributes available
    #    to subfunctions (just a tech detail to overcome a Python 2.x limit)
    class local:
        gph = G
        F = list()
        fxd = xfixed # this is a dictionary of lists not a list of lists
                     # initialize to [] if no fixed vertices
        K = xinit.shape[1]

    # evaluate a random small partial sum_uv (||x_u-x_v||^2-d_uv^2)^2
    #    warning: x is a flat array, needs to be mapped to K x |V|
    def eval_partial_fun(x):
        for v in local.fxd:
            for k in range(local.K):
                x[v*local.K + k] = local.fxd[v][k]
        ret = 0
        for e in local.F:
            ret = ret + (sum((x[e.head.rank*local.K+l] - x[e.tail.rank*local.K+l])**2 for l in range(local.K)) - e.wL**2)**2
        return ret

    # actually do the work
    jac = functools.partial(eval_grad, eval_partial_fun)
    hess = functools.partial(eval_hess, eval_partial_fun)
    convergence = False
    # initial point
    x = xinit.flatten()
    # start the loop
    counter = 0
    print "StochasticNLPSolver(...,", p, ",", close, "): starting"
    eta = 1
    while not convergence:
        counter = counter + 1
        # create a random set of edges
        local.F = list()
        for e in G.E:
            if np.random.uniform(0,1) < p:
                local.F.append(e)
        # we want at least one edge
        if len(local.F) == 0:
            local.F.append(G.E[int(round(np.random.uniform(0,len(G.E))))])
        # save the current point
        y = x
        # compute the next point
        #x = y - np.dot(np.linalg.inv(hess(y)),jac(y))
        x = y - eta * jac(y)
        rc = relativeCloseness(x,y)
        #rc = np.linalg.norm(np.subtract(x,y))
        eta = 0.9*eta
        print "StochasticNLPSolver: closeness =", rc, "at itn", counter
        if rc[0] < close and rc[1] < close:
            convergence = True
    xstar = np.reshape(x, (local.gph.V.card, local.K), order='C')
    return xstar

# MILP model for the DGP in norm 1 (Claudia and Leo)
def Norm1ModelAndSolve(G):
    n = G.V.card
    K = NumberOfDimensions
    t0 = time.time()

    # create problem
    P = pic.Problem()
    # 'big M' bound
    U = sum(e.wL for e in G.E)
    # create variables
    x = P.add_variable('x', (n,K))
    sp = P.add_variable('sp', (n,n))
    sm = P.add_variable('sm', (n,n))
    tp = []
    tm = []
    z  = []
    for k in range(K):
        tp.append(P.add_variable('tp[%d]'%k, (n,n)))
        tm.append(P.add_variable('tm[%d]'%k, (n,n)))
        z.append(P.add_variable('z[%d]'%k, (n,n), vtype = 'integer'))
    # fix non-edge variables
    for i in range(n):
        for j in range(n):
            if not G.isEdge(i,j):
                P.add_constraint(sp[i,j] == 0.0)
                P.add_constraint(sm[i,j] == 0.0)
                for k in range(K):
                    P.add_constraint(tp[k][i,j] == 0.0)
                    P.add_constraint(tm[k][i,j] == 0.0)
                    P.add_constraint(z[k][i,j] == 0.0)
    # bound constraints
    for i in range(n):
        for k in range(K):
            P.add_constraint(x[i,k] > -U)
            P.add_constraint(x[i,k] < U)
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        P.add_constraint(sp[i,j] > 0.0)
        P.add_constraint(sp[i,j] < e.wL)
        P.add_constraint(sm[i,j] > 0.0)
        P.add_constraint(sm[i,j] < e.wL)
        for k in range(K):
            P.add_constraint(tp[k][i,j] > 0.0)
            P.add_constraint(tp[k][i,j] < e.wL)
            P.add_constraint(tm[k][i,j] > 0.0)
            P.add_constraint(tm[k][i,j] < e.wL)
            P.add_constraint(z[k][i,j] > 0.0)
            P.add_constraint(z[k][i,j] < 1)
    # objective function
    P.set_objective('min', sum(sp[i,j] + sm[i,j] for i in range(n) for j in range(n) if G.isEdge(i,j)))
    # constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        # distance
        if e.interval == 1:
            P.add_constraint(e.wL <= sum(tp[k][i,j] + tm[k][i,j] for k in range(K)) - (sp[i,j] - sm[i,j]) <= e.wU)
        else:
            P.add_constraint(sum(tp[k][i,j] + tm[k][i,j] for k in range(K)) - (sp[i,j] - sm[i,j]) == e.wL)
        for k in range(K):
            # relation between t and x
            P.add_constraint(tp[k][i,j] - tm[k][i,j] == x[i,k] - x[j,k])
            # linear complementarity
            P.add_constraint(tp[k][i,j] < U*z[k][i,j])
            P.add_constraint(tm[k][i,j] < U*(1 - z[k][i,j]))
    # zero centroid
    for k in range(K):
        P.add_constraint(sum(x[i,k] for i in range(n)) == 0.0)
    # symmetry
    for i in range(n):
        for j in range(n):
            if i != j:
                P.add_constraint(sp[i,j] == sp[j,i])
                P.add_constraint(sm[i,j] == sm[j,i])
                for k in range(K):
                    P.add_constraint(tp[k][i,j] == tp[k][j,i])
                    P.add_constraint(tm[k][i,j] == tm[k][j,i])
                    P.add_constraint(z[k][i,j] == z[k][j,i])
                    
    t1 = time.time()
    #P.write_to_file('dgp1.lp')
    sol = P.solve(solver = 'cplex', cplex_params = {'mip.display':2})
    xval = np.array(x.value)
    t2 = time.time()
    return (xval, sol)

# MILP model for the DGP in norm 1 to be used in
#   VNS-like heuristic based on limited branching
#   it's a CRAP heuristic that doesn't work and should be ignored (Leo160706)
def Norm1ModelAndSolveLB(G, xstar, radius, timeLimit):
    n = G.V.card
    K = NumberOfDimensions
    t0 = time.time()

    # create problem
    P = pic.Problem()
    # 'big M' bound
    U = sum(e.wL for e in G.E)
    # create variables
    x = P.add_variable('x', (n,K))
    lb = P.add_variable('lb', (n,K), vtype = 'integer')
    sp = P.add_variable('sp', (n,n))
    sm = P.add_variable('sm', (n,n))
    tp = []
    tm = []
    z  = []
    for k in range(K):
        tp.append(P.add_variable('tp[%d]'%k, (n,n)))
        tm.append(P.add_variable('tm[%d]'%k, (n,n)))
        z.append(P.add_variable('z[%d]'%k, (n,n), vtype = 'integer'))
    # fix non-edge variables
    for i in range(n):
        for j in range(n):
            if not G.isEdge(i,j):
                P.add_constraint(sp[i,j] == 0.0)
                P.add_constraint(sm[i,j] == 0.0)
                for k in range(K):
                    P.add_constraint(tp[k][i,j] == 0.0)
                    P.add_constraint(tm[k][i,j] == 0.0)
                    P.add_constraint(z[k][i,j] == 0.0)
    # bound constraints
    for i in range(n):
        for k in range(K):
            P.add_constraint(x[i,k] > -U)
            P.add_constraint(x[i,k] < U)
            P.add_constraint(lb[i,k] > 0)
            P.add_constraint(lb[i,k] < 1)
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        P.add_constraint(sp[i,j] > 0.0)
        P.add_constraint(sp[i,j] < e.wL)
        P.add_constraint(sm[i,j] > 0.0)
        P.add_constraint(sm[i,j] < e.wL)
        for k in range(K):
            P.add_constraint(tp[k][i,j] > 0.0)
            P.add_constraint(tp[k][i,j] < e.wL)
            P.add_constraint(tm[k][i,j] > 0.0)
            P.add_constraint(tm[k][i,j] < e.wL)
            P.add_constraint(z[k][i,j] > 0.0)
            P.add_constraint(z[k][i,j] < 1)
    # objective function
    P.set_objective('min', sum(sp[i,j] + sm[i,j] for i in range(n) for j in range(n) if G.isEdge(i,j)))
    # constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        # distance
        if e.interval == 1:
            P.add_constraint(e.wL <= sum(tp[k][i,j] + tm[k][i,j] for k in range(K)) - (sp[i,j] - sm[i,j]) <= e.wU)            
        else:    
            P.add_constraint(sum(tp[k][i,j] + tm[k][i,j] for k in range(K)) - (sp[i,j] - sm[i,j]) == e.wL)
        for k in range(K):
            # relation between t and x
            P.add_constraint(tp[k][i,j] - tm[k][i,j] == x[i,k] - x[j,k])
            # linear complementarity
            P.add_constraint(tp[k][i,j] < U*z[k][i,j])
            P.add_constraint(tm[k][i,j] < U*(1 - z[k][i,j]))
    # zero centroid
    for k in range(K):
        P.add_constraint(sum(x[i,k] for i in range(n)) == 0.0)
    # symmetry
    for i in range(n):
        for j in range(n):
            if i != j:
                P.add_constraint(sp[i,j] == sp[j,i])
                P.add_constraint(sm[i,j] == sm[j,i])
                for k in range(K):
                    P.add_constraint(tp[k][i,j] == tp[k][j,i])
                    P.add_constraint(tm[k][i,j] == tm[k][j,i])
                    P.add_constraint(z[k][i,j] == z[k][j,i])
    # local branching type constraint
    P.add_constraint(sum(lb[i,k] for i in range(n) for k in range(K)) <= radius)
    # relation between x and lb -- either keeping x[i,k] fixed or not
    for i in range(n):
        for k in range(K):
            P.add_constraint(x[i,k] < xstar[i,k] + U*lb[i,k])
            P.add_constraint(x[i,k] > xstar[i,k] - U*lb[i,k])
    #P.write_to_file('vns1.lp')
    t1 = time.time()
    sol = P.solve(solver = 'cplex', timelimit = timeLimit, cplex_params = {'mip.display':0})
    xval = np.array(x.value)
    t2 = time.time()
    return (xval, sol)

## this method is a VNS-inspired heuristic for solving the DGP in norm 1
##   it's a CRAP heuristic that doesn't work and should be ignored (Leo160706)
def LimitedBranchingVNS(G, maxTime = 120):
    # vertices and dimensions
    n = G.V.card
    K = NumberOfDimensions
    # parameters
    maxUnitSteps = n*K+1  # number of coordinates (max radius)
    stepWidth = 10        # internal loop step
    stepTimeLim = min(n,maxTime/2.0) # max time for subprob (guesswork choice)
    termination = False   # termination flag
    # optimum
    xstar = list()        # best so far
    fstar = myInf         # value of best so far
    fstar1 = fstar        # to decide when to move center to a new solution:
    moveCenter = 0.2      #   when a new sol is found with moveCenter*fstar1
    # the VNS-like heuristic
    t0 = time.time()
    # CPU-based termination
    while not termination:
        # find new center
        center = Matousek(G, K)
        ldecenter = LDE(G, center, 1)
        if ldecenter < fstar:
            fstar = ldecenter
            xstar = center
            print "pydg:LimitedBranchingVNS: Matousek() found new sol with LDE={0:f}".format(fstar)
            if fstar1 == myInf:
                fstar1 = fstar
        for radius in range(stepWidth,maxUnitSteps,stepWidth):
            # try some local search within a limited branching radius 
            (x,err) = Norm1ModelAndSolveLB(G, center, radius, stepTimeLim)
            ldex = LDE(G, x, 1)
            t1 = time.time()
            print "pydg:LimitedBranchingVNS: radius={0:d} CPU={1:f}".format(radius, t1-t0)
            if ldex < fstar:
                fstar = ldex
                xstar = x
                print "pydg:LimitedBranchingVNS: CPLEX found new sol with LDE={0:f}".format(fstar)
                if fstar < moveCenter * fstar1:
                    # change center even though this neighb not fully explored
                    fstar1 = fstar   # re-initialize fstar1 for next time
                    if fstar < myZero:
                        # optimum found
                        termination = True
                    break
            if t1-t0 > maxTime:
                termination = True
                break
    return (xstar, fstar, time.time() - t0)

####### "LIFTED" FORMULATION-BASED METHODS (MATRIX AS VECTOR) ########

## iterative lifted formulation DDP
def IterativeLiftedDDP(G, iterations = 5):
    n = G.V.card
    K = NumberOfDimensions
    nK = n*K
    U = np.identity(nK+1)    
    for i in range(iterations):
        (X,Y) = LiftedDDPModelAndSolve(G,U)
        U = factor(Y)
        lde = LiftedSolutionLDE(G, X)
        print "IterativeLiftedDDP(itn", i+1, "/", iterations, "): rank sol =", np.linalg.matrix_rank(X), "LDE =", lde
        if lde < myZero:
            break
    return X
    

## solve a DGP by lifted SDP, run Barvinok's naive algorithm, then locnlpsolve
def LiftedBarvinok(G, method = "sdp", ddp_itn = 5, locnlp="lbfgs", sdpFormulation=1):
    n = G.V.card
    K = NumberOfDimensions
    nK = n*K
    if method == "sdp":
        sol = LiftedSDPModelAndSolve(G, sdpFormulation)
    elif method == "ddp":
        sol = IterativeLiftedDDP(G, ddp_itn)
    else:
        print "LiftedBarvinok: method should be in {\"sdp\", \"ddp\"}"
        exit
    # clean up solutions from slightly negative eigenvalues and factor it
    X = factor(sol)
    # Barvinok's naive algorithm:
    #    sample from nK-multivariate normal random vector
    y = np.random.multivariate_normal(np.zeros(nK), np.identity(nK))
    #y = np.random.uniform(-1,1,nK)
    #    x = y X
    x1 = np.dot(y,X)
    # re-shape as a realization should be (K x n)
    x2 = np.reshape(x1, (G.V.card,K))
    x = LocalNLPSolver(G,x2,{}, locnlp)
    print "LiftedBarvinok: after Barvinok's alg, LDE =", LDE(G,x2)
    print "LiftedBarvinok: after locnlp", locnlp, ", LDE =", LDE(G,x)    
    return x

## model and solve the U-DDP inner approximation of the lifted SDP using PICOS
def LiftedDDPModelAndSolve(G, U):
    n = G.V.card
    K = NumberOfDimensions
    nK = n*K
    P = pic.Problem()
    # decision variables
    S = P.add_variable('S', (n,n), vtype='symmetric')
    X = P.add_variable('X', (nK,nK), vtype='symmetric')
    Y = P.add_variable('Y', (nK+1,nK+1), vtype='symmetric')
    T = P.add_variable('T', (nK+1,nK+1), vtype='symmetric')
    x = P.add_variable('x', nK)
    # objective functions
    all1 = cvx.matrix(np.ones((n,n)))
    P.set_objective('min', all1 | S)
    ## constraints
    # zero slacks
    for i in range(n):
        for j in range(n):
            if not G.isEdge(i,j):
                P.add_constraint(S[i,j] == 0)                
    # distance constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        if e.interval == 1:
            P.add_constraint(e.wL**2 - S[i,j] <= sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] - 2*X[i*K+l,j*K+l] for l in range(K)) <= e.wU**2 + S[i,j])            
        else:
            P.add_constraint(sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] - 2*X[i*K+l,j*K+l] for l in range(K)) - e.wL**2 <= S[i,j])
            P.add_constraint(-S[i,j] <= sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] - 2*X[i*K+l,j*K+l] for l in range(K)) - e.wL**2)
    # Y is diagonally dominant
    for i in range(nK+1):
        P.add_constraint(sum(T[i,j] for j in range(nK+1) if j != i) <= Y[i,i])
    # Schur complement = U Y U for some DD Y
    Up = cvx.matrix(np.matrix(U))
    UpT = cvx.matrix(np.matrix(U.T))
    P.add_constraint( (( 1 & x.T ) // ( x & X )) == UpT * Y * Up )
    P.add_constraint( Y <=  T )
    P.add_constraint( Y >= -T )
    P.add_constraint( T >= 0 )
    P.add_constraint( S >= 0 )
    # solve and get solution
    P.solve()
    Xval = np.array(X.value)
    Yval = np.array(UpT * Y.value * Up)
    print "LiftedDDPModelAndSolve: LDE =", LiftedSolutionLDE(G, Xval)
    return (Xval, Yval)

## model and solve the lifted SDP formulation using PICOS
def LiftedSDPModelAndSolve(G):
    n = G.V.card
    K = NumberOfDimensions
    nK = n*K

    P = pic.Problem()
    # decision variables
    X = P.add_variable('X', (nK,nK), vtype='symmetric')
    x = P.add_variable('x', nK)
    # objective functions
    ## min sum_{i<j in V} sum_l (X_ilil + X_jljl - 2X_iljl)
    #IJL= [(i,j,l) for i in range(n) for j in range(n) for l in range(K) if i<j]
    ## min sum_{ij in E} sum_l (X_ilil + X_jljl - 2X_iljl)
    #IJL = [(e.tail.rank,e.head.rank,l) for e in G.E for l in range(K)]
    #P.set_objective('min',pic.sum([X[i*K+l,i*K+l] + X[j*K+l,j*K+l] - 2*X[i*K+l,j*K+l] for (i,j,l) in IJL], ('i','j','l'), 'IJL'))
    P.set_objective('find', 1)
    # constraints
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        if e.interval == 1:
            P.add_constraint(e.wL**2 <= sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] -
                                            2*X[i*K+l,j*K+l] for l in range(K)) <= e.wU**2)
        else:
            P.add_constraint(sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] -
                                 2*X[i*K+l,j*K+l] for l in range(K)) == e.wL**2)
    P.add_constraint( (( 1 & x.T ) // ( x & X )) >> 0 )
    # solve and get solution
    P.solve()
    sol = np.array(X.value)
    print "LiftedSDPModelAndSolve: LDE =", LiftedSolutionLDE(G, sol)
    return sol

## verify feasibility of lifted SDP/DDP (matrix) solution
def LiftedSolutionLDE(G,X):
    nK = X.shape[0]
    n = G.V.card
    K = nK / n
    # check solution
    lde = 0
    for e in G.E:
        i = e.tail.rank
        j = e.head.rank
        dij = math.sqrt(sum(X[i*K+l,i*K+l] + X[j*K+l,j*K+l] - 2*X[i*K+l,j*K+l] for l in range(K)))
        if interval == 0:
            newlde = abs(dij - e.wL) / e.wL
        else:
            newlde = max(0, e.wL - dij) + max(0, dij - e.wU)
        if lde < newlde:
            lde = newlde
    return lde

###################### EXTERNAL METHODS ########################

## interface with a dgsol executable in the path -- it might fail
##   if the dgsol output changes format even slightly
def DGSol(G):
    # write .dgsol distance file
    n = G.V.card
    K = NumberOfDimensions
    x = np.zeros((n,K))
    lde = myInf
    mde = myInf
    cpu = 0
    if K != 3:
        print "isomap.py:DGSol(): K must be 3"
        return (x, mde, lde, cpu)
    dgdatf = "isomap_instance.dgsol"
    dgsolf = "dg.sol"
    f = open(dgdatf, "w")
    s = str()
    for e in G.E:
        s = s + e.tail.name + " " + e.head.name + " " + str(e.wL) + " " + str(e.wU) + "\n"  ## based on the fact that e.wL=e.wU for DGP with equations
    f.write(s)
    f.close()
    cmd = "dgsol " + dgdatf
    t0 = time.time()
    msg = os.popen(cmd).read()
    cpu = time.time() - t0
    os.remove(dgdatf)
    x = realizationFromFile(G, dgsolf)
    lde = LDE(G,x)
    mde = MDE(G,x)
    return (x, mde, lde, cpu)


###################### SERVICE ROUTINES (DEBUG) ########################

# plot a realization as points
def plotrlz3D(x):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(x[:,0], x[:,1], x[:,2], c='r', marker='o')
    for i in range(n):
         ax.plot([0, x[i,0]], [0, x[i,1]], [0, x[i,2]], c='b')
    plt.show()
