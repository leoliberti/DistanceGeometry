#!/usr/local/bin/MathematicaScript -script

Get["ApproxRealize.m"];
Get["ProjectionTools.m"];
Get["Projection.m"];

(*
G = Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3, 2 <-> 4, 3 <-> 4}, 
       VertexLabels -> "Name", EdgeWeight -> {2, 2, 1, 1, 1}];
*)

(*
n = 20;
m = Round[N[0.5*n (n-1)/2]];
G = RandomGraph[{n,m}, EdgeWeight -> RandomReal[{0,5}, m]];
A = Normal[WeightedAdjacencyMatrix[G]];
*)

G = SpiralGraph[0.01, 0.6]; 
(* G = SpiralGraph[];*)
xG = VertexCoordinates /. AbsoluteOptions[G, VertexCoordinates];
A = Normal[WeightedAdjacencyMatrix[G]];

(* // test classic MDS on incomplete distance matrix *)
x = ClassicMultiDimensionalScaling[G];
Print["******** classic MDS on partial adj matrix of G"]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]
Print["** project MDS solution to 3D"]
x3 = RandomProjection[x, 3];
Print["Dim Aff x = ", Last[Dimensions[x3]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x3]] ]]

(* // test classic MDS on shortest path matrix  *)
AD = GraphDistanceMatrix[G];
x = ClassicMultiDimensionalScaling[AD];
Print["******** classic MDS on shortest path matrix of G"]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]

(* // test Isomap on weighted graph of local distances  *)
x = Isomap[G, 3];
Print["******** Isomap on G with K=3"]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]

(* // test LLE on weighted graph of local distances *)
x = LocallyLinearEmbedding[G, 3];
Print["******** LLE on G with K=3"]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]

(* // test SPE 
x = StochasticProximityEmbedding[G, 3, 500];
Print["******** SPE on G"]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]
*)

(* // use SPE as a starting point for approx MDS 
B = EuclideanDistanceMatrix[x];
AB = CompletePartialDistanceMatrixBy[A,B];
x = ClassicMultiDimensionalScaling[AB] // Chop;
Print["******** eigenapprox applied to SPE starting point "]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]
*)

(* // refine Isomap with SPE 
x = StochasticProximityEmbedding[G, Isomap[G,3]];
Print["******** refine Isomap realization with SPE "]
Print["Dim Aff x = ", Last[Dimensions[x]]]
Print["relative average error = ", 
      N[RelativePartialEDMError[A, EuclideanDistanceMatrix[x]] ]]
*)

(* // see how SPE's solution can be improved by local descent 
x = StochasticProximityEmbedding[G, 2, 1000];
Print["******** Local minimization from a SPE (1000 itn) starting point"]
Print["SPE error = ", PartialEDMError[A, EuclideanDistanceMatrix[x]]];
Map[(y = DGPSystemApproxLocal[G, x, #]; 
     Print[#, " absolute average error = ",  
	   PartialEDMError[A, EuclideanDistanceMatrix[y]]])&, 
    {"Newton", "QuasiNewton", "LevenbergMarquardt", 
     "ConjugateGradient", "PrincipalAxis"}]
*)

(* // test global optimization methods 
Print["******** Test global optimization methods"]
Map[(y = DGPSystemApproxGlobal[G, 2, #]; 
     Print[#, " absolute average error = ",  
	   PartialEDMError[A, EuclideanDistanceMatrix[y]]])&, 
    {"NelderMead", "DifferentialEvolution", "SimulatedAnnealing", 
     "RandomSearch"}]

Print["******** Test multistart with some local opt methods"]
Map[(y = DGPSystemApproxMultiStart[G, 2, 10, #]; 
     Print[#, " absolute average error = ",  
	   PartialEDMError[A, EuclideanDistanceMatrix[y]]])&, 
    {"Newton", "QuasiNewton", "LevenbergMarquardt", 
     "ConjugateGradient", "PrincipalAxis"}]

y = DGPSystemApproxMultiStart[G, 2];
Print["relative error = ", 
      RelativePartialEDMError[A, EuclideanDistanceMatrix[y]]]
*)
