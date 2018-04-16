#!/usr/local/bin/MathematicaScript -script 

Get["DGPSystem.m"];
Get["CreateRandom.m"];

(* work on the plane *)
K = 2; n = 5; lambda = 0.6;

(* the graph *)
(* G = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 1 <-> 4, 4 <-> 2}, 
       VertexLabels -> "Name", EdgeWeight -> {1, 2, 2, 1, 1}]; *)
G = CreateRandomYES[K,n,lambda];

(* the system *)
Print["DGP System: ", DGPSystem[G,K]];

(* squared system with arbitrarily named variables *)
myVars = Array[x, {Length[VertexList[G]], K}];
S = DGPSystemSquaredVars[G, K, myVars];
Print["DGP System Squared: ", S];

(* fix first three coordinates in realization to remove congruences *)
IncongruentSolutions = {x[1,1]->0, x[1,2]->0, x[2,1]->1};
SI = S /. IncongruentSolutions;

(* solve system symbolically *)
sol = Solve[Apply[And, SI], Flatten[myVars], Reals];
Print["Exact solution: ", sol];

(* solve system numerically *)
Print["Numerical solution 1: ", sol // N];
nsol = NSolve[Apply[And, SI], Flatten[myVars], Reals];
Print["Numerical solution 2: ", nsol];

(* get solutions as list of realizations rather than "rules" *)
sol2 = myVars /. IncongruentSolutions /. sol;

(* plot points in the plane *)
ListPlot[sol2, Joined -> True]

(* display all incongruent solutions *)
Map[GraphPlot[G, VertexLabeling -> True, VertexCoordinateRules -> #] &, sol2];
