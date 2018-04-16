(* write a DGP system *)

DGPSystem[G_Graph, K_Integer, acc_Integer:5] := 
   Block[{x, y, theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 x = Array[y, {Length[vl], K}];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]],acc];
	 theSystem = Map[(
             Sqrt[Sum[ (x[[#[[1]],k]] - x[[#[[2]],k]])^2, {k,K} ]] == 
	     A[[ #[[1]], #[[2]] ]])&, EdgeList[G]
			];
	 Return[theSystem]
	]

DGPSystemSquared[G_Graph, K_Integer, acc_Integer:5] := 
   Block[{x, y, theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 x = Array[y, {Length[vl], K}];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]], acc];
	 theSystem = Map[(
             Sum[ (x[[#[[1]],k]] - x[[#[[2]],k]])^2, {k,K} ] == 
	     (A[[ #[[1]], #[[2]] ]])^2)&, EdgeList[G]
			];
	 Return[theSystem]
	]

DGPSystemAbs[G_Graph, K_Integer, acc_Integer:5] := 
   Block[{x, y, theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 x = Array[y, {Length[vl], K}];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]], acc];
	 theSystem = Map[(
	     EuclideanDistance[ x[[#[[1]]]], x[[#[[2]]]] ] == 
	     A[[ #[[1]], #[[2]] ]])&, EdgeList[G]
			];
	 Return[theSystem]
	]

DGPSystemSquaredAbs[G_Graph, K_Integer, acc_Integer:5] := 
   Block[{x, y, theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 x = Array[y, {Length[vl], K}];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]], acc];
	 theSystem = Map[(
	     (EuclideanDistance[ x[[#[[1]]]], x[[#[[2]]]] ])^2 == 
	     (A[[ #[[1]], #[[2]] ]])^2 )&, EdgeList[G]
			];
	 Return[theSystem]
	]

(* --------- *)

DGPSystemVars[G_Graph, K_Integer, y_, acc_Integer:5] := 
   Block[{theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]], acc];
	 theSystem = Map[(
             Sqrt[Sum[ (y[[#[[1]],k]] - y[[#[[2]],k]])^2, {k,K} ]] == 
	     A[[ #[[1]], #[[2]] ]])&, EdgeList[G]
			];
	 Return[theSystem]
	]

DGPSystemSquaredVars[G_Graph, K_Integer, y_, acc_Integer:5] := 
   Block[{theSystem, A, vl, ip}, 
	 vl = VertexList[G];
	 ip = InversePermutation[vl];
	 A = SetAccuracy[WeightedAdjacencyMatrix[G][[ip,ip]], acc];
	 theSystem = Map[(
             Sum[ (y[[#[[1]],k]] - y[[#[[2]],k]])^2, {k,K} ] == 
	     (A[[ #[[1]], #[[2]] ]])^2)&, EdgeList[G]
			];
	 Return[theSystem]
	]




