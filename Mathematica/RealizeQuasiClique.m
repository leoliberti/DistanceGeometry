(* RealizeQuasiClique.m *)

Get["ApproxRealizeTools.m"];
Get["CayleyMenger.m"];

RealizeQuasiClique[G_Graph?WeightedGraphQ] :=
    Block[{A = Normal[WeightedAdjacencyMatrix[G]], 
	   n, ewl, ewl2, ne, CG, ACG, sol, X1, X, x, w, z, zvars, K},
	  n = Length[A];
	  K = n-2;
	  ewl = EdgeWeightList[A];
	  ne = Count[ewl, x_ /; x>0];
	  If[ne+1 == n (n-1)/2, 
	     zvars = Array[w, {n,n}];
	     z = CayleyMengerVarsSymbolic[A,zvars];
	     sol = Select[Flatten[
		 z /. Solve[CayleyMengerDeterminantSymbolic[A,zvars]==0, z] ], 
			  (#>=0)&];
	     ewl2 = Map[(ReplacePart[ewl, Position[ewl, 0]->#])&, sol];
	     CG = Map[
		 (CompleteGraph[n, EdgeWeight-> 
				ReplacePart[ewl, Position[ewl, 0]->#]])&, sol];
	     X1 = Map[RealizeClique[#]&, CG];
	     X = Map[Take[#, K]&, X1, {2}], 
	     X = -1];
	  Return[X]
	 ]
