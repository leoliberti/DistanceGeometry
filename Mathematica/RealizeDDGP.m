(* RealizeDDGP.m *)

Get["ApproxRealizeTools.m"];
Get["RealizeClique.m"];

On[Assert];

(* pick a random subset of L of cardinality n *)
RandomSet[L_List, n_Integer] :=
    Block[{T = L, S = {}, m = Length[L], x, p},
	  Assert[m >= n];
	  While[Length[T] > m - n, 
		x = (RandomChoice[T, 1])[[1]];
		AppendTo[S, x];
		p = Flatten[Position[T, x]];
		T = Delete[T, p];
	       ];
	  Return[S];
	 ]

(* creates a random DDGP graph (DDGP order 1,...,n) with n>K vertices, 
   0 <= m <= nK-K(K+1)/2 pruning edges and random edge weights, from a 
   random realization with components in [-u,u] *)
DDGPGraph[n_Integer?Positive, K_Integer?Positive, m_Integer, u_Real:5.0] :=
    Block[{x, T, T2, A, B, p, Ui, H, G},
	  Assert[n > K];
	  Assert[m <= n (n-1)/2 - (K (K-1)/2 + K (n-K))];
	  x = RandomReal[{-u,u}, {n, K}];
	  A = EuclideanDistanceMatrix[x];
	  B = ConstantArray[Infinity, {n,n}];
	  Scan[(B[[First[#],Last[#]]] = A[[First[#],Last[#]]];
		B[[Last[#],First[#]]] = A[[Last[#],First[#]]])&, 
	      Flatten[Table[{i,j}, {i,K-1}, {j,Drop[Range[K],i]}], 1]];
	  T = Flatten[Table[{i,j}, {i,1,n-1}, 
			    {j, If[i>=K, Drop[Range[K+1,n],i-K],
				   Range[K+1,n]]}], 1];
	  (* minimal DDGP graph *)
	  For[i = K + 1, i <= n, i++, 
	      Ui = RandomSet[Range[i-1], K];
	      Scan[(B[[#,i]] = A[[#,i]]; B[[i,#]] = A[[i,#]])&, Ui];
	      p = Flatten[Map[(Position[T, #])&, DeleteDuplicates[
		  Flatten[Map[({{i,#}, {#,i}})&, Ui], 1]]], 1];
	      T = Delete[T, p];
	     ];
	  (* add m random pruning edges *)
	  T2 = RandomSet[T, m];
	  Scan[(B[[First[#],Last[#]]] = A[[First[#],Last[#]]];
		B[[Last[#],First[#]]] = A[[Last[#],First[#]]])&, T2];
	  H = WeightedAdjacencyGraph[B];
	  If[K==2||K==3, 
	     G = Graph[Range[n], EdgeList[H], Options[H], 
		       VertexLabels-> "Name", VertexCoordinates->x], 
	     G = Graph[Range[n], EdgeList[H], Options[H], 
		       VertexLabels-> "Name"]
	    ];
	  Return[G];
	 ]

(* creates a random DMDGP graph (DMDGP order 1,...,n) with n>K vertices, 
   0 <= m <= nK-K(K+1)/2 pruning edges and random edge weights, from a 
   random realization with components in [-u,u] *)
DMDGPGraph[n_Integer?Positive, K_Integer?Positive, m_Integer, u_Real:5.0] :=
    Block[{x, T, T2, A, B, p, Ui, H, G},
	  Assert[n > K];
	  Assert[m <= n (n-1)/2 - (K (K-1)/2 + K (n-K))];
	  x = RandomReal[{-u,u}, {n, K}];
	  A = EuclideanDistanceMatrix[x];
	  B = ConstantArray[Infinity, {n,n}];
	  Scan[(B[[First[#],Last[#]]] = A[[First[#],Last[#]]];
		B[[Last[#],First[#]]] = A[[Last[#],First[#]]])&, 
	      Flatten[Table[{i,j}, {i,K-1}, {j,Drop[Range[K],i]}], 1]];
	  T = Flatten[Table[{i,j}, {i,1,n-1}, 
			    {j, If[i>=K, Drop[Range[K+1,n],i-K],
				   Range[K+1,n]]}], 1];
	  (* minimal DMDGP graph *)
	  For[i = K + 1, i <= n, i++, 
	      Ui = Range[i-K,i-1]; (* only difference with DDGPGraph[] *)
	      Scan[(B[[#,i]] = A[[#,i]]; B[[i,#]] = A[[i,#]])&, Ui];
	      p = Flatten[Map[(Position[T, #])&, DeleteDuplicates[
		  Flatten[Map[({{i,#}, {#,i}})&, Ui], 1]]], 1];
	      T = Delete[T, p];
	     ];
	  (* add m random pruning edges *)
	  T2 = RandomSet[T, m];
	  Scan[(B[[First[#],Last[#]]] = A[[First[#],Last[#]]];
		B[[Last[#],First[#]]] = A[[Last[#],First[#]]])&, T2];
	  H = WeightedAdjacencyGraph[B];
	  If[K==2||K==3, 
	     G = Graph[Range[n], EdgeList[H], Options[H], 
		       VertexLabels-> "Name", VertexCoordinates->x], 
	     G = Graph[Range[n], EdgeList[H], Options[H], 
		       VertexLabels-> "Name"]
	    ];
	  Return[G];
	 ]

(* branch-and-prune --- works on DDGP and DMDGP alike *)
RealizeDDGP[G_Graph?WeightedGraphQ, K_Integer?Positive, eps_Real:0.0001] :=
    Block[{X = {}, A = Normal[WeightedAdjacencyMatrix[G]], H, x},
	  H = WeightedAdjacencyGraph[A[[1;;K,1;;K]]];
	  x = Map[PadRight[#,K]&, RealizeClique[H]];
	  X = RealizeDDGPRecursive[A, K, x, K+1, X, eps];
	  Return[X];
	 ]

RealizeDDGPRecursive[A_List?SquareMatrixQ, K_Integer?Positive, 
		     x_List?MatrixQ, i_Integer?Positive, X_List, 
		     eps_Real] :=
    Block[{n = Length[A], Ni, Ui, Upr, W, dst, S, fD, s, y, Y, t},

	  (* adjacent predecessors of i *)
	  Ni = Select[Flatten[Position[A[[i]], z_ /; z > 0]], (# < i)&];
	  Assert[Length[Ni] >= K];
	  (* K closest adjacenct predecessors of i *)
	  Ui = Take[Ni, -K];
	  (* realization of Ui and distances from i to Ui *)
	  W = Map[(x[[#]])&, Ui];
	  dst = Map[(A[[#,i]])&, Ui];
	  (* up to two approximate positions of vertex i *)
	  S = Chop[NextVertexPair[W, dst, K, eps], eps]; 
	  (* j in Upr <==> {j,i} is a pruning edge *)
	  Upr = Take[Ni, Length[Ni]-K];
	 
	  (* prune *)
	  If[Length[Upr] > 0, 
	     s = 1;
	     While[s <= Length[S], 
		   (* compute largest discrepancy with pruning distance *)
		   fD = Max[Map[(Abs[EuclideanDistance[ S[[s]], x[[#]] ] -
				     A[[#,i]] ])&, Upr]];
		   (* if it exceeds eps, prune s *)
		   If[fD > eps, S = Drop[S, {s}], s++];
		  ], 
	     Null];
	 
	  (* branch *)
	  Y = X;
	  For[s = 1, s <= Length[S], s++,
	      (* extend realization x to vertex i *)
	      y = Append[x, S[[s]] ];
	      If[i == n, 
		 (* last vertex: store realization y into pool *)
		 Y = Append[Y,y], 
		 (* not last: recurse *)
		 Y = RealizeDDGPRecursive[A,K,y,i+1,Y,eps] ];
	     ];
	  Return[Y];
	 ]
    
DDGPShow[G_Graph?WeightedGraphQ, K_Integer?Positive, p_Integer:2, X_List:{}] :=
    Block[{Y = X}, 
	  If[X == {}, Y = RealizeDDGP[G,K], Null];
	  If[K == 2, 
	     GraphicsGrid[Partition[Map[(GraphPlot[G, VertexLabeling -> True, 
						   VertexCoordinateRules -> #, 
						   Axes -> False])&,Y], p],
			 Frame -> False],
	     If[K == 3, 
		GraphicsGrid[Partition[
		    Map[(GraphPlot3D[G, VertexLabeling -> True, 
				     VertexCoordinateRules -> #, 
				     Axes -> False])&,Y], p],
			    Frame->False], Null]]
	 ]




(*** oblivion

	  Print["*** placing ", i];
	  Print["  * x ", x];
	  Print["  * adjpred ", Ui];
	  Print["  * x[adjpred] ", W];
	  Print["  * A[adjpred,i] ", dst];
	  Print["  * x_i at ", S];	  
	  Print["  * d(S1,Ui) ", Map[EuclideanDistance[x[[#]],S[[1]]]&, Ui]];
	  Print["  * d(S2,Ui) ", Map[EuclideanDistance[x[[#]],S[[2]]]&, Ui]];
	  If[Apply[And, Map[(Element[#, Reals]) &, Flatten[S]]], Null, Quit[]];

	  Print["  * Upr ", Upr];
	  Print["  * d(S1,Upr) ", Map[EuclideanDistance[x[[#]],S[[1]]]&, Upr]];
	  Print["  * d(S2,Upr) ", Map[EuclideanDistance[x[[#]],S[[2]]]&, Upr]];

*)