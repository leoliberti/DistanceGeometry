(*********** Cayley-Menger determinant **********)

(* test to establish whether a given matrix is square *)
SquareMatrixQ[A_List] := MatrixQ[A] && 
    First[Dimensions[A]] == Last[Dimensions[A]];

(********** numeric versions **********)

(* compute the CM matrix from a (complete) distance matrix *)
CayleyMengerMatrix[A_List?SquareMatrixQ] := 
    Prepend[Map[(Prepend[#, 1]) &, Map[(#^2)&, A]],
	    Prepend[ConstantArray[1, Length[A]], 0]];

(* compute the CM determinant from a (complete) distance matrix *)
CayleyMengerDeterminant[A_List?SquareMatrixQ] := 
    Det[CayleyMengerMatrix[A]];

(* compute the CM volume from a (complete) distance matrix *)
CayleyMengerVolume[A_List?SquareMatrixQ, K_Integer?Positive] := 
    Sqrt[(-1)^(K)/(2^(K+1) ((K+1)!)^2) CayleyMengerDeterminant[A]];

(* (complete) graph equivalents *)
CayleyMengerMatrix[G_Graph?WeightedGraphQ] := 
    CayleyMengerMatrix[Normal[WeightedAdjacencyMatrix[G]]];

CayleyMengerDeterminant[G_Graph?WeightedGraphQ] := 
    Det[CayleyMengerMatrix[G]];

CayleyMengerVolume[G_Graph?WeightedGraphQ, K_Integer?Positive] := 
    Sqrt[(-1)^(K)/(2^(K+1) ((K+1)!)^2) CayleyMengerDeterminant[G]];

(************** symbolic versions ***********)

(* get the index pairs of the missing distances *)
CayleyMengerSymbolic[A_List?SquareMatrixQ, var_List?MatrixQ] :=
    Cases[Position[A, 0], {i_, j_} /; i < j];

(* get the list of variables used in the symbolic CM matrix *)
CayleyMengerVarsSymbolic[A_List?SquareMatrixQ, var_List?MatrixQ] :=
    Map[var[[First[#], Last[#]]]&, CayleyMengerSymbolic[A,var]];

(* make a CM matrix from an incomplete distance matrix, using var symbols
   to replace missing distances *)
CayleyMengerMatrixSymbolic[A_List?SquareMatrixQ, var_List?MatrixQ] :=
    Block[{B = A}, Scan[(B[[Last[#], First[#]]] = 
	   B[[First[#], Last[#]]] = var[[First[#], Last[#]]]) &,
	  CayleyMengerSymbolic[A,var]];
     CayleyMengerMatrix[B]];

CayleyMengerDeterminantSymbolic[A_List?SquareMatrixQ, var_List?MatrixQ] :=
    Det[CayleyMengerMatrixSymbolic[A, var]];

CayleyMengerVolumeSymbolic[A_List?SquareMatrixQ, K_Integer?Positive, 
			   var_List?MatrixQ] := 
    Sqrt[(-1)^(K)/(2^(K+1) ((K+1)!)^2) CayleyMengerDeterminantSymbolic[A, var]];

(* graph versions *)
CayleyMengerMatrixSymbolic[G_Graph?WeightedGraphQ, var_List?MatrixQ] :=
    CayleyMengerMatrixSymbol[Normal[WeightedAdjacencyMatrix[G]], var];

CayleyMengerDeterminantSymbolic[G_Graph?WeightedGraphQ, var_List?MatrixQ] :=
    Det[CayleyMengerMatrixSymbolic[G, var]];

CayleyMengerVolumeSymbolic[G_Graph?WeightedGraphQ, K_Integer?Positive, 
			   var_List?MatrixQ] := 
    Sqrt[(-1)^(K)/(2^(K+1) ((K+1)!)^2) CayleyMengerDeterminantSymbolic[G, var]];

(**************** unit simplex **************)

UnitMatrix[K_Integer?Positive] := 
    Table[1, {i, 1, K}, {j, 1, K}];

SimplexDistanceMatrix[K_Integer?Positive] := 
    UnitMatrix[K] - IdentityMatrix[K];
