#!/usr/local/bin/MathematicaScript -script 

Get["NextVertex.m"];

K = 3;

x = {{0,0}, {0,1}, {1,0}};
d = {Sqrt[2], 1, 1};
Print["x_4 = ", NextVertexUnique[x,d,K]]

x = {{0,0,0}, {0,1,0}, {1,0,0}};
d = {1, 1, 1};
Print["x_4 in ", NextVertexPair[x,d,K]]

