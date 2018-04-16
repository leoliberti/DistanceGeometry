(* Mathematica lists matrices by row - always use transpose *)

EdgeOrientation[g_Graph] := 
    Graph[Map[{#[[1]],#[[2]]}, EdgeList[G]], DirectedEdges -> True];

RigidityMatrix[g_Graph, x_List?MatrixQ] := 
    Block[{ K = Dimensions[x][[2]], n = Dimensions[x][[1]], 
	    EL = EdgeList[G], m, R, ed, secant }, 

	  m = Length[EL];
	  R = ConstantArray[0, {m,n K}];

	  For[i = 1, i <= m, i++, 
	      ed = EL[[i]];
	      For[j = 1, j <= n, j++, 
		  If[j == ed[[1]],
		     secant = x[[j]] - x[[ed[[2]]]],
		     If[j == ed[[2]],
			secant = x[[j]] - x[[ed[[1]]]],
			secant = ConstantArray[0,K]
		       ]
		    ];
		  For[k = 1, k <= K, k++,
		      R[[i,K (j-1) + k]] = secant[[k]]
		     ]
		 ]
	     ];
	  Return[R];
	 ]
    