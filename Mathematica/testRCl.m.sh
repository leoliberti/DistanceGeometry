#!/usr/local/bin/MathematicaScript -script 

Get["RealizeClique.m"];

G = CompleteGraph[4, EdgeWeight -> {1, 1, 1, 1, 1, 1}, 
  VertexLabeling -> True, VertexLabels -> "Name"];
x = RealizeClique[G];
GraphPlot3D[G, VertexLabeling -> True, VertexCoordinateRules -> x, 
	    Axes -> True, AxesStyle -> Gray];
Print[x]


n = 10;
G = CompleteGraph[n, EdgeWeight -> ConstantArray[1, n (n-1)/2]];
x = RealizeClique[G];
xProj = Map[(Take[#, 3]) &, x];
GraphPlot3D[G, VertexLabeling -> True, 
 VertexCoordinateRules -> xProj, Axes -> True, AxesStyle -> Gray]
Print[x]