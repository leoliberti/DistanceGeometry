#!/usr/local/bin/MathematicaScript -script 

Get["RealizeDDGP.m"];
Get["RealizeTrilaterative.m"];
Get["VertexOrder.m"];
Get["ProjectionTools.m"];

(*
G = Graph[{1<->2, 1<->3, 1<->5, 1<->6, 1<->8, 2<->3, 
	2<->4, 2<->5, 2<->6, 3<->4, 3<->5, 3<->8, 
	4<->5, 4<->7, 5<->6, 5<->7, 6<->7, 6<->8, 7<->8}];
alpha = TrilaterationOrder[G, 2];
Print[alpha];
*)

G = SpiralGraph[0.4, 2.2];
alpha2 = Timing[TrilaterationOrder[G, 2]];
Print[alpha2]
alpha3 = Timing[TrilaterationOrder[G, 3]];
Print[alpha3]
alpha4 = Timing[TrilaterationOrder[G, 4]];
Print[alpha4]
alpha5 = Timing[TrilaterationOrder[G, 5]];
Print[alpha5]

(*
G = DDGPGraph[20, 3, 5];
pi = RandomPermutation[20];
vtxReplRules = Map[(# -> Permute[VertexList[G], pi][[#]]) &, VertexList[G]];
--- THIS WON'T WORK!! produces graphs with missing vertices in their list! ---
H = VertexReplace[G, vtxReplRules];
Timing[TrilaterationOrder[G, 3]];
*)

(*
G = DDGPGraph[6, 2, 1];
omega = ContiguousTrilaterationOrder[G, 2];
Print[EdgeList[G]]
Print[omega]
*)

