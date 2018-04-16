(* ApproxRealizeTools.m  - various tools used by ApproxRealize.m *)

(* the DGP function (sum of squared errors from DGPSystem) *)
DGPFunction[G_Graph?WeightedGraphQ, x_List?MatrixQ] :=
    Total[Map[(SquaredEuclideanDistance[x[[First[#]]],x[[Last[#]]]] - 
	       Normal[WeightedAdjacencyMatrix[G]][[First[#],Last[#]]])^2&, 
	      EdgeList[G]]]

(* form the distance matrix for a realization *)
EuclideanDistanceMatrix[x_List?MatrixQ] :=
    Block[{n = First[Dimensions[x]]},
	  Partition[Map[(EuclideanDistance[x[[First[#]]], x[[Last[#]]]]) &, 
			Flatten[Table[{i,j}, {i,n}, {j,n}], 1]], n]]

(* test to establish whether a given matrix is square *)
SquareMatrixQ[A_List] := MatrixQ[A] && 
    First[Dimensions[A]] == Last[Dimensions[A]]

(* Edge list weights from weighted adjacency matrix *)
EdgeWeightList[A_List?SquareMatrixQ] :=
    Map[(A[[First[#],Last[#]]])&, 
	Flatten[Table[{i,j}, {i,Length[A]-1}, {j,i+1,Length[A]}],1]]

(* mean absolute error between a partial matrix A and a matrix B, computed
   on the nonzero components of A *)
PartialEDMError[A_List?SquareMatrixQ, B_List?SquareMatrixQ] :=
    Mean[Map[(Abs[A[[First[#], Last[#]]] - B[[First[#], Last[#]]]]) &, 
	     Position[A, x_ /; x > 0]]]

(* mean relative error between a partial matrix A and a matrix B, computed
   on the nonzero components of A *)
RelativePartialEDMError[A_List?SquareMatrixQ, B_List?SquareMatrixQ] :=
    Mean[Map[((Abs[A[[First[#], Last[#]]] - B[[First[#], Last[#]]]]) /
	      A[[First[#], Last[#]]]) &, Position[A, x_ /; x > 0]]]

(* mean of A_ij / B_ij computed on nonzero A_ij's *)
PartialEDMRatio[A_List?SquareMatrixQ, B_List?SquareMatrixQ] :=
    Mean[Map[(B[[First[#], Last[#]]] / A[[First[#], Last[#]]])&, 
	     Position[A,x_ /; x>0]]]

(* complete a partial distance matrix by a complete one: nonzero
   values of the partial one are carried over, zeroes are overwritten *)
CompletePartialDistanceMatrixBy[A_List?SquareMatrixQ, B_List?SquareMatrixQ] :=
    Block[{n = First[Dimensions[A]], n2 = First[Dimensions[B]], AB}, 
	  Assert[n == n2];
	  AB = Partition[Map[(If[A[[First[#], Last[#]]] > 0, 
				 A[[First[#], Last[#]]], 
				 B[[First[#], Last[#]]]]) &, 
			     Flatten[Table[{i, j}, {i, n}, {j, n}], 1]], n];
	  Return[AB];
	 ]

(* translate a realization so that x[[1]] = 0 *)
RealizationTranslateToOrigin[x_List?MatrixQ] :=
    Map[(x[[#]] - x[[1]])&, Range[Length[x]]]	  

(* find the centroid of a set of points *)
Centroid[x_List?MatrixQ] := Total[x] / Length[x]

(* find the distance from centroid to farthest point *)
CentroidRadius[x_List?MatrixQ] := 
    Max[Map[EuclideanDistance[Centroid[x], #]&, x]]
