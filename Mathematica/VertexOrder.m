(* VertexOrder.m *)

(* relabel vertices of a graph: 
   V must be the image of a permutation of the graph vertices, 
   ordered as {1,...,n} *)
VertexRelabel[G_Graph, V_List] :=
    VertexReplace[G, Map[(# -> L[[#]]) &, Range[Length[VertexList[G]]]]]

(* find a K-trilateration order (a DDGP order) *)
TrilaterationOrder[G_Graph, K_Integer?Positive] :=
    Block[{V = VertexList[G], alpha = {}, P, CC}, 
	  P = Subsets[V,{K}];
	  While[alpha == {} && Length[P]>0,
		CC = First[P];
		P = Rest[P];
		alpha = TrilaterationOrder[G, CC];
	       ];
	  Return[alpha];
	 ]

(* find a |C|-trilateration order from a given initial clique C *)
TrilaterationOrder[G_Graph, CC_List] :=
    Block[{alpha = CC, K = Length[CC], V = VertexList[G], n, W, a, amaxW, v},
	  n = Length[V];
	  a = ConstantArray[0, n];
	  W = Select[V, (!MemberQ[alpha,#])&];
	  Scan[(a[[#]] = Length[Intersection[AdjacencyList[G, #], alpha]])&, W];
	  While[Length[W]>0, 
		amaxW = Max[Map[a[[#]]&, W]]; 
		v = Intersection[Flatten[Position[a, amaxW]], W][[1]];
		If[a[[v]] < K, alpha = {}; Break[], Null];
		AppendTo[alpha, v];
		Scan[(a[[#]]=a[[#]]+1)&, Intersection[AdjacencyList[G,v],W]];
		W = Drop[W, Flatten[Position[W, v]]];
	       ];
	  Return[alpha];
	 ]

(* find a contiguous K-trilateration order (a DMDGP order) *)
ContiguousTrilaterationOrder[G_Graph, K_Integer?Positive] :=
    Block[{n = Length[VertexList[G]], V = VertexList[G], t, y, z, NR,
	   uniqueRank, uniqueVertex, initialClique, contiguousAdjPred, 
	   constraints, bounds, sol},
	  If[ConnectedGraphQ[G] == True,
	     t = Array[x, {n, n}];
	     NR = Range[n];
	     uniqueRank = Map[(Sum[t[[#, i]], {i, n}] == 1) &, NR];
	     uniqueVertex = Map[(Sum[t[[v, #]], {v, n}] == 1) &, V];
	     initialClique = Map[(Sum[t[[u, j]], 
				      {u, AdjacencyList[G, First[#]]}, 
				      {j, 1, Last[#] - 1}] >= 
				  (Last[#]-1) t[[First[#], Last[#]]]) &, 
				 Flatten[Table[{v,i}, {v,V}, {i,2,K}], 1]];
	     contiguousAdjPred = Map[(Sum[t[[u, j]], 
					  {u, AdjacencyList[G, First[#]]}, 
					  {j, Last[#] - K, Last[#] - 1}] >= 
				      K t[[First[#], Last[#]]]) &, 
				     Flatten[Table[{v,i},{v,V},{i,K+1,n}],1]];
	     constraints = Apply[And, 
				 Union[uniqueRank, uniqueVertex, 
				       initialClique, contiguousAdjPred]];
	     bounds = Apply[And, 
			    Map[(t[[First[#], Last[#]]] >= 0 && 
				 t[[First[#], Last[#]]] <= 1)&, 
				Flatten[Table[{v,i}, {v,V}, {i,n}],1]]];
	     sol = Minimize[{0, constraints && bounds }, Flatten[t], Integers];
	     y = t /. sol[[2]];
	     If[sol[[1]] == Infinity, z = {}, 
		z = Map[(Sum[v y[[v,#]], {v,V}])&, NR]]; ,
	     z = {}];
	  Return[z];
	 ]

