function [ka, smax, smin] = sfs_ocm(smax, smin, invA, Duu, Dvv, H, Nomegain)
% Optimizaion of second derivative magnitudes and signs
%
% -- Input
% smax : Large surface second derivative sign. [Nomegain, 1]
% smin : Small surface second derivative sign. [Nomegain, 1]
% invA : Inverse of matrix A. [Nomega, Nomega]
% Duu : Matrix for differential oparator d^2/du^2. [Nomegain, Nomega]
% Dvv : Matrix for differential oparator d^2/dv^2. [Nomegain, Nomega]
% H : Diagonal matrix with elements (1-aniso). [Nomegain, Nomegain]
% Nomegain : number of pixels in omegain2d
%
% -- Output
% ka : Optimized second derivative magnitudes. [Nomegain, 1]
% smax : Optimized large surface second derivative sign. [Nomegain, 1]
% smin : Optimized small surface second derivative sign. [Nomegain, 1]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

Smax = spdiags(smax,0,Nomegain,Nomegain);
Smin = spdiags(smin,0,Nomegain,Nomegain);
S = Dvv'*Smax + Duu'*H*Smin;
R = (speye(Nomegain)+H*H) - S'*invA*S;

ka = R \ ones(Nomegain,1);
smax = sign(ka).*smax;
smin = sign(ka).*smin;
ka = abs(ka);
ka = ka./mean(ka);

end

