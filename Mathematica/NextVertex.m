(* NextVertex.m -- positioning the next vertex *)
(* x = x_1,...,x_K; d = d_1,...,d_K [where d_i=d_{i,K+1}] *)

NextVertexUnique[x_List?MatrixQ, d_List, K_Integer?Positive, 
		 eps_Real:0.00001] :=
    Block[{MyVars, MySystem, z},
	  MyVars = Array[y, K-1];
	  MySystem = Map[(2 (x[[#]]-x[[K]]).MyVars == 
			  Norm[x[[#]]]^2 - Norm[x[[K]]]^2 -
			  d[[#]]^2 + d[[K]]^2)&, 
			 Range[K-1]];
	  sol = Solve[Apply[And, MySystem], MyVars];
	  z = Flatten[MyVars /. sol];
	  If[Abs[EuclideanDistance[z, x[[K]]]-d[[K]]] > eps, z = {}, Null];
	  Return[z]
	 ]

NextVertexPair[x_List?MatrixQ, d_List, K_Integer?Positive, 
	      eps_Real:0.00001] :=
    Block[{MyVars, MyMatrix, MyRHS, y, sol, S},
	  MyVars = Array[y, K];
	  MyMatrix = Map[(2 (x[[#]]-x[[K]]))&, Range[K-1]];
	  MyRHS = Map[(Norm[x[[#]]]^2 - Norm[x[[K]]]^2 - 
		       d[[#]]^2 + d[[K]]^2)&, Range[K-1]];
	  MySystem = {MyMatrix.MyVars == MyRHS,
		      Sum[(MyVars[[k]] - x[[1,k]])^2, {k,K}] == d[[1]]^2};
	  sol = Solve[Apply[And, MySystem], MyVars];
	  S = Reverse[MyVars /. sol];
	  (* see if it's really just one point *)
	  If[Length[S] == 2 &&
	     PossibleZeroQ[Total[Chop[MapThread[(#1-#2)&, S], eps]]], 
	     S = S[[{1}]], Null];
	  (* if solutions are complex, ditch them *)
	  If[Apply[And, Map[(Element[#, Reals]) &, Flatten[S]]], Null, S = {}];
	  Return[S];
	]

(*
MaxLDE[z_List, x_List, d_List, K_Integer] :=
    Max[Map[(Abs[EuclideanDistance[z, x[[#]]] - d[[#]]])&, Range[K]]]
    
SumLDE[z_List, x_List, d_List, K_Integer] :=
    Total[Map[(Abs[EuclideanDistance[z, x[[#]]] - d[[#]]])&, Range[K]]]
*)

(**************** oblivion *************)

(*
	 x0 = LinearSolve[MyMatrix, MyRHS];
         ker = NullSpace[MyMatrix];
	 kdir = Array[t, Length[ker]];
         v = x0 + kdir.ker;
*)
