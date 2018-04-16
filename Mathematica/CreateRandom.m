(******* create random DGP instances ******)

Get["ProjectionTools.m"];

(* create a random DGP instance with given number of vertices and 
   edge density lambda; instance is NO with probability 1 *)
CreateRandomNO[n_Integer?Positive, lambda_Real] :=
  Block[{G, m}, 
        m = Round[N[lambda*(n*(n-1)/2)]];
	G = RandomGraph[{n,m}, EdgeWeight -> RandomReal[{0,10},m]];
	(* GraphicsRow[{GraphPlot[G, VertexLabeling -> True],
		     MatrixPlot[Normal[WeightedAdjacencyMatrix[G]]]}] *)
	Return[G];
       ]

(* create a YES random DGP instance in given dimension K, with
   given number of vertices and edge density lambda *)
CreateRandomYES[K_Integer?Positive, n_Integer?Positive, lambda_Real] :=
  Block[{G, x, A, S}, 
	x = RandomReal[{-10,10}, {n, K}];
	A = EuclideanDistanceMatrix[x];
	S = Abs[SkewSymmetrize[Floor[RandomReal[{lambda,1+lambda}, {n,n}]]]];
	A = S*A /. 0. -> Infinity;
	G = WeightedAdjacencyGraph[A, VertexCoordinates -> x];
	Return[G];
       ]
