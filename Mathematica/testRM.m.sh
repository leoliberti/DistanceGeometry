#!/usr/local/bin/MathematicaScript -script 

Get["RigidityMatrix.m"];

(* the graph *)
G = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 1 <-> 4, 4 <-> 2}, 
       VertexLabels -> "Name", EdgeWeight -> {1, 2, 2, 0.5, 0.5}];

(* the realization *)
(*x = {{0, 1}, {0, -1}, {-1, 0}, {1, 0}};*)
x = PropertyValue[{G, #}, VertexCoordinates] & /@ VertexList[G]; (*automatic*)
Print["E(G) = ", Map[{#[[1]],#[[2]]}&, EdgeList[G]]];
Print["x_G = ", x];

R = RigidityMatrix[G,x];
Print["Rank of rigidity matrix = ", MatrixRank[R]];
