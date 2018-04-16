(* RealizeClique *)

Get["NextVertex.m"];

RealizeClique[G_Graph?WeightedGraphQ] :=
    Block[{i, j, n, x, S, K},
	  n = Length[VertexList[G]];
	  Assert[n > 1];
	  K = n - 1;
	  x = ConstantArray[0, {n, K}];
	  Scan[(x[[1,#]] = 0)&, Range[K]];
	  Scan[(x[[2,#+1]] = 0)&, Range[K-1]];
	  x[[2,1]] = WeightedAdjacencyMatrix[G][[1,2]];
	  For[i = 3, i <= n, i++, 
	      S = NextVertexPair[Table[Take[x[[j]], i-1], {j, 1, i-1}],
		  Table[WeightedAdjacencyMatrix[G][[j,i]], {j, 1, i-1}], 
		  i-1];
	      If[S == {}, Return[{}], Null];
	      x[[i]] = PadRight[First[S], K];
	     ];
	  Return[x];
	 ]

RealizeClique[A_List?SquareMatrixQ] :=
    RealizeClique[WeightedAdjacencyGraph[A]]

    