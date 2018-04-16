#!/usr/local/bin/MathematicaScript -script

Get["ApproxRealize.m"];
Get["Projection.m"];

(* // Johnson-Lindenstrauss and random projections to small dimensions 
n = 100;
bnd = 10;
x0 = RandomReal[{-bnd,bnd}, {n,n}];
A = EuclideanDistanceMatrix[x0];
Print["******** Johnson-Lindenstrauss"]
Amean = Mean[Select[Flatten[A], Positive]];
x = ClassicMultiDimensionalScaling[A];
xJL = JohnsonLindenstrauss[x, 0.3];
If[IntegerQ[xJL] && xJL < 0, 
   Print["low distortion embedding in higher dim than ambient space", xJL],
   xdim = Last[Dimensions[xJL]];
   AJL = EuclideanDistanceMatrix[xJL];
   Print["Number of points = ", First[Dimensions[A]]]
   Print["JL embedding dimension = ", xdim]
   Print["JL relative per-distance error = ", RelativePartialEDMError[A, AJL]] ]
Print["******** Random projection to R^3"]
xproj = RandomProjection[A, 3];
Aproj = EuclideanDistanceMatrix[xproj];
Print["RP relative per-distance error = ", RelativePartialEDMError[A, Aproj]]
*)

(* // local Nash device *)
Print["******** Local Nash device"]
G = SpiralGraph[];
xG = VertexCoordinates /. AbsoluteOptions[G, VertexCoordinates];
x = LiftRandomRotate[xG, 100];
A = EuclideanDistanceMatrix[x];
Apartial = Normal[WeightedAdjacencyMatrix[G]];
k = Round[Mean[Map[(Total[#])&, AdjacencyMatrix[G]]]];
Print[k]
y = LocalNashDevice[x, k];
Print[Dimensions[y]]
Ay = EuclideanDistanceMatrix[y];
Print["Nash relative error = ", RelativePartialEDMError[Apartial, Ay]]
Print["Nash relative error with full distance matrix = ", 
      RelativePartialEDMError[A, Ay]]
