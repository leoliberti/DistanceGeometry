% read the 'tiny' instance
[P,K] = readAMPLdatpedm('dat/tiny_gph.dat');
% realize it using SDP
[x,ret] = sdprealize(K,P);
% show the realization
showrealization(x);
% realize it using DDP
[y,ret] = ddprealize(K,P);
% realize it using iterative DDP (5 iterations)
[z,ret] = iterddp(K,P,5);
% show all realizations at once
X(1,:,:) = x;
X(2,:,:) = y;
X(3,:,:) = z;
showrealizations(X);
