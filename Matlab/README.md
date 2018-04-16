Matlab code for Distance Geometry
Code by Leo Liberti (c) 2014-2017

Ipartialerror.m          compute avg partial l_2 rlz error, interval version
Ipedmerror.m             compute avg pEDM error, interval version
IreadAMPLdatpedm.m       read an AMPL .dat file out of gph2dat, interval version
Irlzerror.m              compute avg realization error, interval version
Irndddgp.m               generate random iDDGP instances
Irnddmdgp.m              generate random iDMDGP instances
Irndpedm.m               generate random iDGP instances
Isdprealize.m            solve iDGP instances via SDP (needs Yalmip)
Iyajima.m                solve iDGP via SDP using Yajima's formul. (Yalmip)
alignment.m              interface for comparing alignment modulo part refl
alignpartial.m           use absor.m to align initial segments of 2 realiz.
alignprefl.m             align two rlz modulo rot, transl, and part refl
alignrealization.m       use absor.m to align 2 realizations
allcliques.m		         find all K-cliques in the pEDM P (used by ddgporder)
allprefl.m               generate all partial refl of given rlz & pruning group
bpbestlinsys.m           used by bpnext
bpbestquadratic.m        used by bpnext
bplinsys.m               used by bpbestlinsys
bpnext.m                 compute the next point positions (used by branchprune)
bpquadratic.m            used by bpbestquadratic
branchprune.m            branch & prune
chirality.m              compute chirality of realization
countchr.m               number of given characters in string
dgpeq.m                  ? [to be used within fmincon?]
dist2gram.m              transform an EDM into a Gram matrix
disterror.m              =norm(A-B)
ddgporderfromclique.m    DDGP order from starting clique (used by ddgporder)
ddgporder.m              computes DDGP orders 
dmdgporder.m             computes kDMDGP orders (needs CPLEX)
dscrprnmat.m             partitions DDGP into discretization and pruning edges
eps2zero.m               zero all rlz components < given epsilon
eucldist.m               Euclidean distance matrix from realization
floydwarshall.m          compute all shortest paths
gbulinsys.m              geometric build-up algorithm (used by bplinsys)
givens.m                 Givens rotation matrix
givens2.m                2D Givens rotation matrix (used by givens)
graph2pedm.m             turn weighted graph into a partial EDM
hhreflmat.m              Householder reflection matrix
isddgp.m                 check whether the pEDM P is a K DDGP instance
isdmdgp.m                check whether the pEDM P is a K-DMDGP instance
isomap.m                 the Isomap method
ispartialreflection.m    ascertain whether 2 rlz are part refl of each other
jll.m                    Johnson-Lindenstrauss projections
lde.m			               compute max partial rlz error over partial EDM entries
losemptyedm.m            remove -1 columns/rows from pEDM
losemptyrlz.m            remove origin from realizations
mde.m               	   compute avg error between rlz and pEDM
mds.m                    multidimensional scaling
partialerror.m           compute avg partial rlz error over partial EDM entries
partialreflection.m      partial reflection of rlz at vertex v
pca.m                    principal component analysis
pedm2adj.m               adjacency matrix from a pEDM
pedmdraw.m               draw 2D/3D realizations with all graph edges
pedmerror.m              compute avg pEDM error
preflapprox.m            greedy partial reflection alignment between 2 rlz
pruninggroup.m           compute the pruning group from a pEDM
readAMPLdatpedm.m        read an AMPL _gph.dat file out of gph2dat, get pEDM
readAMPLdatrlz.m         read an AMPL _rlz.dat or _sol.dat with a realization
realizeclique.m          compute one realization of an n-clique in R^{n-1}
reflection.m             compute reflection and apply it to rlz
reflectionmatrix.m       reflection matrix 
reflectionoperator.m     compute reflection operator between 2 rlz
rndNMRedm.m              generate random DGP instance with systematic error
rndddgp.m                generate random DDGP instance
rnddmdgp.m               generate random DMDGP instance
rndpedm.m                generate random DGP instance
rotationmatrix.m         rotation matrix
sdprealize.m             solve DGP instance using SDP (needs Yalmip)
sdprealize2.m            solve DGP instance using SDP, alternate form (Yalmip)
showrealization.m        display a 2D/3D realization    
showrealizations.m       display a set of 2D/3D realizations
spe.m                    stochastic proximity embedding
sqeucldist.m             squared Euclidean distance matrix from realization
symmfromupper.m          symmetrize a matrix (copy UR triangle to LL)
transl0centroid.m        translate rlz so centroid is zero
