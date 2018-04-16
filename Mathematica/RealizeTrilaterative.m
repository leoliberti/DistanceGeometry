(* RealizeTrilaterative.m *)

Get["NextVertex.m"];
Get["RealizeClique.m"];

(* creates a trilaterative graph (order 1,...,n) with
   random edge weights, from a reandom realization with components in [-U,U] *)
TrilaterativeGraph[n_Integer?Positive, K_Integer, U_Real:5.0] :=
    Block[{x, el, dst, G},
	  Assert[n > K+1];
	  x = RandomReal[{-U,U}, {n, K}];
	  el = Flatten[Union[Table[{j<->i}, {i,2,K+1}, {j,1,i-1}], 
			     Table[{j<->i}, {i,K+2,n}, {j,i-K-1,i-1}]]];
	  dst = Map[(EuclideanDistance[x[[First[#]]],x[[Last[#]]]])&, el];
	  If[K==2||K==3, 
	     G = Graph[Range[n], el, EdgeWeight->dst, VertexLabels->"Name", 
		       VertexCoordinates->x],
	     G = Graph[Range[n], el, EdgeWeight->dst, VertexLabels->"Name"]
	    ];
	  Return[G];
	 ]

(* assumes the trilaterative order is 1,...,n *)
RealizeTrilaterative[G_Graph?WeightedGraphQ, K_Integer?Positive] :=
    Block[{A = WeightedAdjacencyMatrix[G], n, B, H, x, i, j},
	  n = Length[A];
	  B = A[[1;;K+1,1;;K+1]];
	  H = WeightedAdjacencyGraph[B];
	  x = RealizeClique[H];
	  Scan[AppendTo[x, ConstantArray[0,K]]&, Range[n-K-1]];
	  For[i = K+2, i <= n, i++, 
	      x[[i]] = NextVertexUnique[
		  Table[x[[j]], {j, i-K-1, i-1}],
		  Table[A[[j,i]], {j, i-K-1, i-1}], 
		  K+1];
	      For[j = 1, j < i-K-1, j++, 
		  If[A[[j,i]]>0 && EuclideanDistance[x[[i]],x[[j]]] != A[[j,i]],
		     x[[i]] = {}; Break[], Null];
		 ];
	      If[x[[i]] == {}, x = -1, Null];
	     ];
	  Return[x];
	 ]
