(* Projection.m - projection methods *)

Get["ApproxRealizeTools.m"];
Get["ApproxRealize.m"];
Get["ProjectionTools.m"];

(****** random standard Gaussian projection on a given dimension ******)
(* Warning: matrix.vector scalar product is actually vector.matrix,
   hence the interchange of the small and large dimension, i.e.
   the projection matrix is tall and thin, not large and squat *)
RandomProjection[x_List?MatrixQ, K_Integer?Positive] :=
    Block[{R, H = Last[Dimensions[x]]},
	  R = RandomVariate[NormalDistribution[0, Sqrt[1/K]], {H,K}];
	  Return[ N[x.R] ]
	 ]

(*** J-L method, a la [Indyk, Motwani] (dense projection matrix) 
     the constant C (default=1.8) is from [Venkatasubramanian, Wang] ***)
JohnsonLindenstrauss[x_List?MatrixQ, eps_Real:0.3] :=
    Block[{K, n = First[Dimensions[x]], C = 1.8},
	  K = Round[(C/(eps^2)) Log[n]];
	  Return[RandomProjection[x, K]]
	 ]

(************** Local Nash device [Bartal et al.] ****************)
LocalNashDevice[x_List?MatrixQ, k_Integer, eps_Real:0.3] := 
    Block[{A = EuclideanDistanceMatrix[x], 
	   H, r, sA, sigmaC, omegaC, amplitudeC, dist2borderC, 
	   Cl, nCl, rCl, Clx, lCl, ii, i, j, h, y, nuC,
	   n = First[Dimensions[x]], K = Last[Dimensions[x]] },

	  H = Round[(0.4 /(eps^2)) Log[k]];
	  (* neighbourhood radius for given k *)
	  r = Max[Map[(Last[Take[Sort[sA = A[[#]]], k]])&, Range[n]]]; 
	  (* partition data and find cluster radii *)
	  nCl = Round[N[n/k]];
	  Cl = FindClusters[x -> Range[n], nCl];
	  Clx = Map[(Map[x[[#]]&, Cl[[#]]])&, Range[nCl]];
	  rCl = Map[CentroidRadius[Clx[[#]]]&, Range[nCl]];
	  rCl = Map[(If[Abs[#]<eps^2,r,#])&, rCl];
	  (* compute embedding parameters *)
	  sigmaC = Map[(N[(1 / (rCl[[#]] eps)) Log[k]])&, Range[nCl]];
	  omegaC = Map[(RandomVariate[NormalDistribution[0, 1], {H,K}])&, 
		      Range[nCl]];
	  nuC = Map[(RandomVariate[BernoulliDistribution[1/2], H])&, 
		    Range[nCl]];
	  (* compose embedding *)
	  y = ConstantArray[0, {n, 3 H}];
	  For[j = 1, j <= nCl, j++,
	      (* for each cluster *)
	      lCl = Length[ Cl[[j]] ];
	      For[i = 1, i <= lCl, i++,
		  (* for each point in the cluster *)
		  ii = Cl[[j,i]];
		  dist2borderC = Min[Map[(A[[ ii, #]])&, 
					 Complement[Range[n], Cl[[j]] ]]] ;
		  amplitudeC = Min[{dist2borderC, 1/sigmaC[[j]]}];
		  For[h = 1, h <= H, h++, 
		      (* for odd dimension indices in the projection *)
		      y[[ii, 2 h - 1]] = 
		      amplitudeC Cos[sigmaC[[j]] omegaC[[j,h]].x[[ii]]];
		     ];
		  For[h = 1, h <= H, h++, 
		      (* for even dimension indices in the projection *)
		      y[[ii, 2 h]] = 
		      amplitudeC Sin[sigmaC[[j]] omegaC[[j,h]].x[[ii]]];
		     ];
		  For[h = 1, h <= H, h++,
		      y[[ii, 2 H + h]] = 
		      Sqrt[eps] nuC[[j,h]] dist2borderC;
		     ];
		 ];
	     ];
	  Return[y]
	 ]
