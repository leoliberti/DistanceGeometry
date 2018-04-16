#!/usr/local/bin/MathematicaScript -script 

Get["RealizeComplete.m"];

(* the graph *)
G = Graph[{1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 2 <-> 4, 3 <-> 4}, 
       VertexLabels -> "Name", 
	  EdgeWeight -> {1, 2, 1, 2, 1, (Sqrt[15]+Sqrt[3])/2}];
K = 2;
x0 = {{0,0}, {1,0}, {1/2, Sqrt[15]/2}};
x = RealizeComplete[G, 2, x0];
Print[x];

(* realize a 5-clique in R^2 *)
K5 = CompleteGraph[5, VertexLabels -> "Name", 
		   EdgeWeight -> {1,1,1,1,Sqrt[2],2,Sqrt[2],Sqrt[2],2,Sqrt[2]}];
x0 = {{0,0}, {Sqrt[2]/2, Sqrt[2]/2}, {Sqrt[2]/2, -Sqrt[2]/2}};
x = RealizeComplete[K5, 2, x0];
Print[x];
GraphPlot[K5, VertexLabeling -> True, VertexCoordinateRules -> x]
