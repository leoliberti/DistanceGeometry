#!/usr/local/bin/MathematicaScript -script 

Get["RealizeDDGP.m"];
Get["RealizeTrilaterative.m"];

(*
G2 = DDGPGraph[7, 2, 2];
X = RealizeDDGP[G2, 2];
Print[Length[X]]
*)

G = Graph[{1<->2, 1<->3, 1<->7, 2<->3, 2<->4, 2<->5, 2<->6, 3<->4, 
	3<->5, 3<->6, 5<->7, 6<->7}, 
    EdgeWeight -> {9.53541, 3.11071, 6.88687, 7.34597, 7.58441, 
	7.8516, 1.16467, 6.28141, 7.96418, 8.47708, 3.18685, 5.16973}];
X = RealizeDDGP[G, 2];
Print[Length[X]]
