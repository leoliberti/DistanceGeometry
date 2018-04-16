(* ProjectionTools.m *)

Get["ApproxRealizeTools.m"];

(* copy (- upper triangle) into lower triangle *)
SkewSymmetrize[A_List?SquareMatrixQ] := 
    Block[{n = Length[A], B = A}, 
	  Scan[(B[[Last[#], First[#]]] = -B[[First[#], Last[#]]]) &, 
	       Select[Flatten[Table[{i, j}, {i, n}, {j, n}], 1],
		      (First[#] < Last[#])&]];
	  Return[B]
	 ]

(* see http://en.wikipedia.org/wiki/Rotation_matrix 
   [maths/papers/wikipedia-rotation_matrix.pdf]
   "Skew parameters via Cayley's Formula *)
CayleyRotationMatrix[A_List?SquareMatrixQ] :=
    Block[{Askew = SkewSymmetrize[A], 
	   Idm = IdentityMatrix[Length[A]]}, 
	  Return[(Idm + Askew) . Inverse[Idm - Askew]]
	 ]

(* Lift realization to n dimensions and rotate by a random matrix *)
LiftRandomRotate[x_List?MatrixQ, n_Integer] :=
    Block[{R, liftx, y}, 
	  liftx = Map[(PadRight[#, n])&, x];
	  R = CayleyRotationMatrix[RandomReal[{-1,1}, {n,n}]];
	  y = liftx . R;
	  Return[y]
	 ]

(* produce a weighted graph with local distances *)
SpiralGraph[step_:0.02, localRadius_:0.7] :=
    Block[{x, n, ESet, WSet, G},
	  x = Map[SpiralGraphFunction, Range[3, 7, step]];
	  n = Length[x];
	  ESet = Select[Flatten[Table[{i, j}, {i, n}, {j, n}], 1], 
			(First[#] < Last[#] && 
			 EuclideanDistance[x[[First[#]]], x[[Last[#]]]] < 
			 localRadius)&];
	  WSet = Map[(EuclideanDistance[x[[First[#]]], x[[Last[#]]]])&, ESet];
	  G = Graph[Range[n], ESet, VertexCoordinates -> x, EdgeWeight -> WSet];
	  Return[G]
	 ]
SpiralGraphFunction[t_] := 
    {Sin[10 t] + Sqrt[100 - t^2], Cos[10 t] + Sqrt[100 - t^2], 
     10/t} + RandomReal[{-0.2, 0.2}, 3]

(* return the vertex coordinates from a graph *)
Realization[G_Graph] :=
    VertexCoordinates /. AbsoluteOptions[G, VertexCoordinates]

