# DistanceGeometry
Research codes for distance geometry

This file details the distribution of my research codes in Distance Geometry (DG). Most of my studies in DG concern the DG Problem (DGP), which is the inverse problem to the statement "given some points in R^K, compute some of the pairwise distances".

The formal definition of the DGP is as follows: given an integer K>0 and a weighted undirected simple graph G=(V,E,d), where d assigns nonnegative values to the edges in E, is there a realization x of the vertices V such that, for each edge {u,v} in E, x satisfies ||x_u - x_v|| = d_uv ?
 
Among the best known variants, we have the Euclidean DGP (EDGP), where the norm Euclidean (and the distance function d maps to square distance values), the Molecular DGP (MDGP), which is a EDGP where K=3, and the interval DGP (iDGP) where d_uv is an interval, and the condition to satisfy is ||x_u - x_v|| in d_uv.

The DGP is a decision problem which is known to be strongly NP-hard even for fixed K (and in particular K=1) [Saxe 1978]. It is in NP for K=1 but its NP membership status is unknown for other values of K.

The DGP is an appropriate abstraction of many problems in science and engineering. For K=1, we have the Clock Synchronization Problem, where the absolute times in a network of clocks must be computed given their phase differences; for K=2 we have the Sensor Network Localization Problem, where the position of wireless sensors is computed starting from some of the pairwise distances; for K=3 we have the Protein Conformation Problem, where the geometric structure of a protein is computed from a set of inter-atomic distances measured using Nuclear Magnetic Resonance.

The research codes in this repository are:
- in Mathematica (these accompany the book [Liberti, Lavor, Euclidean Distance Geometry: an Introduction, 2017])
- in MATLAB (these have been developed to test ideas in several research papers)
- in Python (these have been developed to test ideas in several research papers).

