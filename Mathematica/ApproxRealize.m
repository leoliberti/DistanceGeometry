(* ApproxRealize.m -- Approximate realizations of graphs *)

Get["ApproxRealizeTools.m"];

(******* classic MultiDimensional Scaling (MDS) [Cox&Cox] *******)

(* given a complete EDM, find the realization *)
ClassicMultiDimensionalScaling[A_List?SquareMatrixQ, eps_Real:0.000001] := 
    Block[{n = Length[A], J, B, lambda, V, x, pos, l1, l2},
	  (* positive definite eigenspace of -1/2 J A^2 J *)
	  J = IdentityMatrix[n] - ConstantArray[1., {n,n}]/n;
	  B = J.(-A*A/2).J;
	  lambda = N[Eigenvalues[B]] // Chop;
	  pos = Position[lambda, x_ /; x > eps];
	  V = Extract[N[Eigenvectors[B]], pos];
	  x = RealizationTranslateToOrigin[
	      Transpose[Extract[Sqrt[lambda], pos]*V]];
	  Return[x];
	 ]

(* treat a partial matrix as if it were complete and hope for the best *)
ClassicMultiDimensionalScaling[G_Graph?WeightedGraphQ] := 
    Block[{x, l1, l2, A = Normal[WeightedAdjacencyMatrix[G]]},
	  n = Length[A];
	  x = ClassicMultiDimensionalScaling[A];
	  (* distances are just scaled, scale back *)
	  l1 = Map[(EuclideanDistance[x[[First[#]]], x[[Last[#]]]]) &, 
		   Tuples[Range[n], 2]];
	  l2 =  Flatten[A];
	  mm = Mean[Flatten[Map[(l1[[#]]/l2[[#]]) &, 
				Position[l2, x_ /; x > 0]]]];
	  Return[x/mm];
	 ]

(********* Isomap applied to weighted graphs *********)

(*** Isomap [Tenenbaum, de Silva, Langford] on graphs 
     (local distances are weighted edges) *****)
Isomap[G_Graph?WeightedGraphQ, K_Integer?Positive] :=
    ApproxMultiDimensionalScaling[GraphDistanceMatrix[G], K];

(* given a complete EDM, find an approximate realization in K dim *)
ApproxMultiDimensionalScaling[A_List?SquareMatrixQ, K_Integer?Positive] :=
    Block[{n = Length[A], J, B, lambda, pos, y},
	  J = IdentityMatrix[n] - ConstantArray[1., {n,n}]/n;
	  B = J.(-A*A/2).J;
	  lambda = Take[N[Eigenvalues[B]], K] // Chop;
	  V = Take[N[Eigenvectors[B]], K];
	  x = RealizationTranslateToOrigin[Transpose[Sqrt[lambda]*V]];
	  Return[x];
	 ]

(******* Locally Linear Embedding applied to weighted graphs *******)

LocallyLinearEmbedding[G_Graph?WeightedGraphQ, K_Integer?Positive] :=
    Block[{A = Normal[WeightedAdjacencyMatrix[G]], 
	   x, y, n, d, W, Wi, v1, Zi, nZ, Ci, IZ, B, 
	   neighbSize, rankCtrl, pos},
	  (* find high-dimensional realization of G *)
	  x = ClassicMultiDimensionalScaling[GraphDistanceMatrix[G]];
	  n = Length[x];
	  d = Last[Dimensions[x]];
	  (* compute optimal weights *)
	  W = ConstantArray[0, {n,n}];
	  neighbSize = Max[Map[(Count[A[[#]], x_ /; x>0])&, Range[n]]];
	  If[neighbSize > d, rankCtrl = 1, rankCtrl = 0];
	  For[i = 1, i <= n, i++,
	      pos = Flatten[Position[A[[i]], x_ /; x > 0]];
	      nZ = Length[pos];
	      Zi = ConstantArray[0, {nZ, d}];
	      Scan[(Zi[[#]] = x[[ pos[[#]] ]] - x[[i]])&, Range[nZ]];
	      Ci = Zi . Transpose[Zi];
	      Ci += 0.001 Tr[Ci] rankCtrl IdentityMatrix[nZ];
	      Wi = LinearSolve[Ci, ConstantArray[1, nZ]];
	      Wi /= Total[Wi];
	      For[j = 1, j <= Length[pos], j++, 
		  W[[i, pos[[j]] ]] = Wi[[j]];
		 ];
	     ];
	  (* compute projection *)
	  IZ = IdentityMatrix[n] - W;
	  B = IZ. Transpose[IZ];
	  y = Transpose[Reverse[Take[Eigenvectors[B], {n-K, n-1}]]];
	  Return[y]
	 ]

(**************** SPE heuristic ***************)

(* taken from [Agrafiotis et al, in Mucherino et al.] p. 293 (alg 7) *)
StochasticProximityEmbedding[G_Graph?WeightedGraphQ, K_Integer?Positive, 
			     Steps_Integer:1000] :=
    Block[{n = Length[VertexList[G]], 
	   A = Normal[WeightedAdjacencyMatrix[G]], r, x0, x},
	  (* random initial embedding *)
	  r = N[Total[Flatten[A]]/K];
	  x0 = RealizationTranslateToOrigin[RandomReal[{-r,r}, {n, K}]];
	  x = StochasticProximityEmbedding[G, x0, Steps];
	  Return[RealizationTranslateToOrigin[x]];
	 ]

(* SPE used as an improvement on an initial point *)
StochasticProximityEmbedding[G_Graph?WeightedGraphQ, x0_List?MatrixQ, 
			     Steps_Integer:1000] :=
    Block[{n = Length[VertexList[G]], K = Last[Dimensions[x0]],
	   A = Normal[WeightedAdjacencyMatrix[G]], 
	   x, i, j, decr, eps, cycles, lambda = 0.9},
	  (* constants *)
	  cycles = n;
	  decr = N[0.89 / n];
	  eps = N[Mean[Select[Flatten[A], Positive]]/100];
	  x = x0;
	  (* spe *)
	  For[i = 1, i <= cycles, i++,
	      (* run a SPE iteration *)
	      For[j = 1, j <= Steps, j++, 
		  x = SPEIteration[G, A, x, K, eps, lambda];
		 ];
	      (* adjust learning rate *)
	      lambda -= decr;
	     ];
	  Return[x];
	 ]


SPEIteration[G_Graph?WeightedGraphQ, A_List?MatrixQ, x_List?MatrixQ,
	     K_Integer?Positive, eps_Real, lambda_Real] :=
    Block[{t, u, v, dst, rst, delta, y, zero = 0.000000001},
	  (* spe *)
	  t = RandomChoice[EdgeList[G]];
	  u = First[t];
	  v = Last[t];
	  dst = EuclideanDistance[x[[u]], x[[v]]];
	  rst = A[[u,v]];
	  y = x;
	  delta = If[Abs[dst - rst] > eps,
		     ((lambda/2) (rst-dst) / (dst+zero)) (y[[u]] - y[[v]]), 
		     ConstantArray[0, K]];
	  y[[u]] += delta;
	  y[[v]] -= delta;
	  Return[y];
	 ]

(******* minimize squared error: local and global opt. methods **********)

DGPSystemApproxLocal[G_Graph?WeightedGraphQ, 
		     StartingPoint_List?MatrixQ, 
		     MethodName_String:"Newton"] :=
    Block[{x, y, z, sol, n, K}, 
	  n = First[Dimensions[StartingPoint]];
	  K = Last[Dimensions[StartingPoint]];
	  y = Array[z, {n, K}];
	  sol = 
	  FindMinimum[DGPFunction[G,y], 
		      Map[({y[[First[#], Last[#]]], 
			    StartingPoint[[First[#], Last[#]]]}) &, 
			  Flatten[Table[{i, k}, {i, n}, {k, K}], 1]],
		      Method -> MethodName
		     ];
	  x = RealizationTranslateToOrigin[y /. Last[sol]];
	  Return[x]
	 ]

DGPSystemApproxGlobal[G_Graph?WeightedGraphQ, K_Integer?Positive,
		      MethodName_String:"DifferentialEvolution"] :=
    Block[{x, y, z, sol, n}, 
	  n = Length[VertexList[G]];
	  y = Array[z, {n, K}];
	  sol = 
	  NMinimize[{DGPFunction[G,y]}, 
		    Map[(y[[First[#], Last[#]]]) &, 
			Flatten[Table[{i, k}, {i, n}, {k, K}], 1]],
		      Method -> MethodName
		     ];
	  x = RealizationTranslateToOrigin[y /. Last[sol]];
	  Return[x]
	 ]

DGPSystemApproxMultiStart[G_Graph?WeightedGraphQ, K_Integer?Positive, 
			 IterationLimit_Integer:10, 
			  MethodName_String:"Newton"] :=
    Block[{x, y, z, xstar, f, fstar, r, i, 
	   A = WeightedAdjacencyMatrix[G], 
	   n = Length[VertexList[G]]},
	  y = Array[z, {n, K}];
	  r = N[Total[Flatten[A]]/K];
	  xstar = RandomReal[{-r,r}, {n, K}];
	  fstar = DGPFunction[G, xstar];
	  For[i = 1, i <= IterationLimit, i++,
	      x = RandomReal[{-r,r}, {n, K}];
	      sol = FindMinimum[DGPFunction[G, y],
				Map[({y[[First[#],Last[#]]],
				      x[[First[#],Last[#]]]})&,
				    Flatten[Table[{i, k}, {i, n}, {k, K}], 1]],
			  Method -> MethodName
			 ];
	      f = First[sol];
	      x = RealizationTranslateToOrigin[y /. Last[sol]];
	      If[f < fstar, xstar = x; fstar = f, Null];
	     ];
	  Return[x];
	 ]
    
