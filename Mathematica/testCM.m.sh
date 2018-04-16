#!/usr/local/bin/MathematicaScript -script 

Get["CayleyMenger.m"]

G = Graph[{1 <-> 2, 1 <-> 3, 2 <-> 3, 2 <-> 4, 3 <-> 4}, 
	  VertexLabels -> "Name", EdgeWeight -> {2, 2, 1, 1, 1}]
A = Normal[WeightedAdjacencyMatrix[G]]
Print[A]

Print[CayleyMengerMatrix[A]]
Print[CayleyMengerDeterminant[A]]
Print[CayleyMengerVolume[A,2]]

y = Array[x, Dimensions[A]]
sol = CayleyMengerVarsSymbolic[A, y] /. 
 Solve[CayleyMengerDeterminantSymbolic[A, y] == 0, 
  CayleyMengerVarsSymbolic[A, y]]

Print[sol]
Print[N[sol]]
