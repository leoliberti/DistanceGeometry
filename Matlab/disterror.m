% e = disterror(A,B)
%   compute difference ||A-B||

function e = disterror(A,B)
  e = norm(A-B)
end