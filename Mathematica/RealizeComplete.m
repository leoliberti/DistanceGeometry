(* RealizeComplete *)

Get["NextVertex.m"];

RealizeComplete[G_Graph?WeightedGraphQ, K_Integer?Positive, x0_List?MatrixQ] :=
    Block[{i, j, n, x, A = WeightedAdjacencyMatrix[G]},
	  Assert[Length[x] == K+1];
	  n = Length[VertexList[G]];
	  x = ConstantArray[0, {n, K}];
	  Scan[(x[[#]] = x0[[#]])&, Range[K+1]];
	  For[i = K+2, i <= n, i++, 
	      x[[i]] = NextVertexUnique[
		  Table[x[[j]], {j, i-K-1, i-1}],
		  Table[A[[j,i]], {j, i-K-1, i-1}], 
		  K+1];
	      For[j = 1, j < i-K-1, j++, 
		  If[EuclideanDistance[x[[i]],x[[j]]] != A[[j,i]],
		     x[[i]] = {}; Break[], Null];
		 ];
	      If[x[[i]] == {}, x = -1, Null];
	     ];
	  Return[x];
	 ]
