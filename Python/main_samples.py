############# main from formulations.py #############

if len(sys.argv) < 2:
    print "error, need .gph file on cmd line"
    exit(1)

# set number of dimensions
NumberOfDimensions = 3

# read graph from data file
inputFile = sys.argv[1]
G = Graph(edgesFromFile(inputFile))
filename = os.path.basename(inputFile)

# SDP formulation
sdpFormulation = 1
if len(sys.argv) >= 3:
    sdpFormulation = int(sys.argv[2])

if sdpFormulation < 1 or sdpFormulation > 3:
    print "error, sdpFormulation must be in {1,2,3}"
    print "    sdpFormulation = 1: [Push-and-pull]"
    print "    sdpFormulation = 2: [Ye 2003]"
    print "    sdpFormulation = 3: [Barvinok 1995]"
    exit(2)
    
x = dict()
cpu = dict()
allmethods = ["rndrlz", "ms", "isomap", "isomap+nlp", "sdp", "sdp+nlp", "sdp+barvinok", "ddp", "ddp+barvinok", "iterddp+barvinok", "ddp+nlp", "liftsdp", "liftddp", "iterddp"]

localNLPsolver = "lbfgs"
#localNLPsolver = "ipopt"

#mymethods = ["sdp+barvinok"]
mymethods = ["ddp+barvinok"]
#mymethods = ["iterddp+barvinok"]

for m in mymethods:
    if m == "rndrlz":
        t0 = time.time()
        x[m] = RandomRealization(G, NumberOfDimensions, 100) #100=bound
        cpu[m] = time.time() - t0
    if m == "ms":
        t0 = time.time()
        x[m] = MultiStart(G, 20, method=localNLPsolver) #20=iterations
        cpu[m] = time.time() - t0
    if m == "isomap":
        t0 = time.time()
        x[m] = IsoMap(G)
        cpu[m] = time.time() - t0
    if m == "isomap+nlp":
        t0 = time.time()
        x[m] = LocalNLPSolver(G,IsoMap(G),{}, method=localNLPsolver)
        cpu[m] = time.time() - t0
    if m == "sdp":
        t0 = time.time()
        x[m] = SDP(G,sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "sdp+nlp":
        t0 = time.time()
        x[m] = LocalNLPSolver(G, SDP(G,sdpFormulation), {}, method=localNLPsolver)
        cpu[m] = time.time() - t0        
    if m == "sdp+barvinok":
        t0 = time.time()
        x[m] = Barvinok(G, "sdp", locnlp=localNLPsolver, sdpFormulation=sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "ddp":
        t0 = time.time()
        x[m] = DDLP(G, objtype=sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "iterddp":
        t0 = time.time()
        x[m] = PCA(IterativeDDLP(G, DDPItn, objtype=sdpFormulation), NumberOfDimensions)
        cpu[m] = time.time() - t0
    if m == "ddp+nlp":
        t0 = time.time()
        x[m] = LocalNLPSolver(G, DDLP(G, objtype=sdpFormulation), {}, method=localNLPsolver)
        cpu[m] = time.time() - t0
    if m == "ddp+barvinok":
        t0 = time.time()
        x[m] = Barvinok(G, "ddp", ddp_itn=DDPItn, locnlp=localNLPsolver, sdpFormulation=sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "iterddp+barvinok":
        t0 = time.time()
        x[m] = Barvinok(G, "iterddp", ddp_itn=DDPItn, locnlp=localNLPsolver, sdpFormulation=sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "liftsdp":
        t0 = time.time()
        x[m] = LiftedBarvinok(G, "sdp", locnlp=localNLPsolver, sdpFormulation=sdpFormulation)
        cpu[m] = time.time() - t0
    if m == "liftddp":
        t0 = time.time()
        x[m] = LiftedBarvinok(G, "ddp", ddp_itn=DDPItn, locnlp=localNLPsolver, sdpFormulation=sdpFormulation)
        cpu[m] = time.time() - t0

# out: one method per line
print "instance method lde mde cpu"
for m in mymethods:
    print filename, m, LDE(G, x[m]), MDE(G, x[m]), cpu[m]

# out: one line
print "OUT", filename, G.V.card, len(G.E), 
for m in mymethods:
    print "{0:0.2f}".format(LDE(G, x[m])),
for m in mymethods:
    print "{0:0.2f}".format(MDE(G, x[m])),
for m in mymethods:
    print "{0:0.2f}".format(cpu[m]),
print

#G.plot(xlift)

################# main from dg_isomap.py #################
if len(sys.argv) < 2:
    print "error, need .gph file on cmd line"
    exit(1)

# set number of dimensions
NumberOfDimensions = 3

# read graph from data file
inputFile = sys.argv[1]
G = Graph(edgesFromFile(inputFile))
filename = os.path.basename(inputFile)

x = dict()
cpu = dict()
allmethods = ["rndrlz", "ms", "isomap", "isomap+nlp", "isomapheur", "sdp", "sdp+nlp", "rndproj", "spe", "barvinoksdp"]
alllocnlp = ["ipopt", "lbfgs", "cg", "newtoncg", "powell", "neldermead"]

## best overall for both speed and quality: ipopt
#localNLPsolver = "ipopt"
#localNLPsolver = "lbfgs"
#localNLPsolver = "newtoncg"
#localNLPsolver = "cg"
#localNLPsolver = "powell"
#localNLPsolver = "neldermead"
localNLPsolver = "dgsol"
#localNLPsolver = "leastsq"

if len(sys.argv) >= 3:
    if sys.argv[2] in alllocnlp:
        localNLPsolver = sys.argv[2]
    
mymethods = ["barvinoksdp", "sdp"]

for m in mymethods:
    if m == "rndrlz":
        t0 = time.time()
        x[m] = RandomRealization(G, NumberOfDimensions, 100) #100=bound
        cpu[m] = time.time() - t0
    if m == "spe":
        t0 = time.time()
        x[m] = StochasticProximityEmbedding(G, RandomRealization(G, NumberOfDimensions, 100), p=0.001) #100=bound
        cpu[m] = time.time() - t0
    if m == "ms":
        t0 = time.time()
        x[m] = MultiStart(G, 20, method=localNLPsolver) #20=iterations
        cpu[m] = time.time() - t0
    if m == "rndproj":
        t0 = time.time()
        x[m] = RandomProjection(G)
        cpu[m] = time.time() - t0
    if m == "isomap":
        t0 = time.time()
        x[m] = IsoMap(G)
        cpu[m] = time.time() - t0
    if m == "isomap+nlp":
        t0 = time.time()
        x[m] = LocalNLPSolver(G,IsoMap(G),{}, method=localNLPsolver, maxitn=nlpMaxIter-2)
        cpu[m] = time.time() - t0
    if m == "isomap+spe":
        t0 = time.time()
        x[m] = StochasticProximityEmbedding(G,IsoMap(G), p=0.01)
        cpu[m] = time.time() - t0
    if m == "isomapheur":
        t0 = time.time()
        x[m] = IsoMapHeuristic(G, locnlp=localNLPsolver)
        cpu[m] = time.time() - t0
    if m == "sdp":
        t0 = time.time()
        x[m] = SDP(G)
        cpu[m] = time.time() - t0
    if m == "sdp+nlp":
        t0 = time.time()
        x[m] = LocalNLPSolver(G, SDP(G), {}, method=localNLPsolver, maxitn=nlpMaxIter)
        cpu[m] = time.time() - t0        
    if m == "barvinoksdp":
        t0 = time.time()
        bvkbnd = int(round((math.sqrt(8*len(G.E) + 1)-1)/2))
        (X,rk) = BarvinokSDP(G,bvkbnd)
        print "BarvinokSDP: bound for dimension =", bvkbnd
        print "BarvinokSDP: solution rank =", rk
        y = PCA(X,NumberOfDimensions)
        x[m] = LocalNLPSolver(G,y,{}, method=localNLPsolver, maxitn=nlpMaxIter-2)
        #x[m] = y
        cpu[m] = time.time() - t0

# out: one method per line
print "instance method lde mde cpu"
for m in mymethods:
    print filename, m, LDE(G, x[m]), MDE(G, x[m]), cpu[m]

# out: one line
print "OUT", localNLPsolver, filename, G.V.card, len(G.E), 
for m in mymethods:
    print LDE(G, x[m]),
for m in mymethods:
    print MDE(G, x[m]),
for m in mymethods:
    print cpu[m],
print

#G.plot(xlift)


################# main from bp.py ##################

# set number of dimensions
NumberOfDimensions = 3

# read graph from data file
inputFile = sys.argv[1]
dmdgp = instanceFromFile(inputFile)
G = dmdgp.graph

x = RealizeClique(G, 3)
print x
print distanceMatrix(x)

## test a few methods
#xrr = RandomRealization(G, NumberOfDimensions, 10)
#xms = MultiStart(G, 20)
#xiso = IsoMap(G)
#xipiso = IPOpt(G,IsoMap(G),{})
#xsdp = SDP(G)
#xipsdp = IPOpt(G,SDP(G),{})
#xddlp = DDLP(G)
#xipddlp = IPOpt(G,DDLP(G),{})

#print "###", 
#print LDE(G,xrr)
#print LDE(G,xms)
#print LDE(G,xiso)
#print LDE(G,xipiso)
#print LDE(G,xddlp)
#print LDE(G,xipddlp)
#print LDE(G,xsdp)
#print LDE(G,xipsdp)

G.plot(x)

## orders
#(alpha,yesno,nat2top,top2nat) = DDGPOrder(G)
#if not yesno:
#    print "NO"
#else:
#    print alpha
#    (theta,phi) = bondAnglesCosines(G,alpha)
#    print theta
#    print phi

#print bondAnglesCosines(G, range(G.V.card))

## next1D
#p = [4]                    
#Z = next1D(p,1.5)

# next2D
#p = [np.array([0,0]), np.array([1,0])]
#Z = next2D(p, [1,1])
#print p, Z

## next3D
#p = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
#dists = [1.0,1.0,1.0]
#p = [np.array([1.3,-0.2,0.6]), np.array([0.1,-1.7,-0.3]), np.array([1,-1.2,1])]
#dists = [2.0,2.0,2.0]
#Z = next3D(p,dists)
#Z = next3Dtrig(p,dists)
#Z = next(p,dists,3)
# print Z
# print "{0:d} {1:d} {2:f}".format(0,1,np.linalg.norm(np.subtract(p[0],p[1])))
# print "{0:d} {1:d} {2:f}".format(0,2,np.linalg.norm(np.subtract(p[0],p[2])))
# print "{0:d} {1:d} {2:f}".format(1,2,np.linalg.norm(np.subtract(p[1],p[2])))
# for z in Z:
#     print "---"
#     for i in range(3):
#         print "{0:d} 3 {1:f}".format(i,np.linalg.norm(np.subtract(z,p[i])))
