#!/usr/local/bin/MathematicaScript -script 

Get["RealizeTrilaterative.m"];

G = TrilaterativeGraph[7,2];
x = RealizeTrilaterative[G,2];
Print[x]

G2 = TrilaterativeGraph[7,3];
x = RealizeTrilaterative[G2,3];
Print[x]

