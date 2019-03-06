function [ka, A2] = sfs_gd(ka, smax, smin, A2, B, Duu, Dvv, Dall, H, Q, Nomegain, parm_gd)
% Gradient descent for optimizing second derivative magnitudes
%
% -- Input
% ka : Curvature magnitudes. [Nomegain, 1]
% smax : Large surface second derivative sign. [Nomegain, 1]
% smin : Small surface second derivative sign. [Nomegain, 1]
% A2 : A2 = Duu'*P^2*Duu + Dvv'*P^2*Dvv + 2*Duv'*P^2*Duv. [Nomega, Nomega]
% B : Matrix for boundary condition. [Nomega, Nomega]
% Q : Matries used in gradient descent. {Nomegain, 1}
% Nomegain : number of pixels in omegain2d
%
% -- Output
% ka : Optimized second derivative magnitudes. [Nomegain, 1]
% smax : Optimized large surface second derivative sign. [Nomegain, 1]
% smin : Optimized small surface second derivative sign. [Nomegain, 1]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

p = 1./ka; % inverse of second derivative magnitude. p = 1./ka, P = inv(Ka)
Smax = spdiags(smax,0,Nomegain,Nomegain);
Smin = spdiags(smin,0,Nomegain,Nomegain);
S = Dvv'*Smax + Duu'*H*Smin;
Ntype_dxdy = size(Dall.Duu0,1);

for repeat = 1 : parm_gd.N
    invASp = A2\(S*p);
    tmp = S'*invASp;
    ka_new = zeros(Nomegain,1);
    for i = 1 : Nomegain
        ka_new(i) = ka(i) + parm_gd.gamma * p(i)^2 * ( p(i)*(invASp'*Q{i}*invASp) - tmp(i) );
    end
    P = spdiags(1./ka_new,0,Nomegain,Nomegain);
    p = full(diag(P));
    ka = 1./p;
    A2 = B;
    for itype = 1 : Ntype_dxdy
        A2 = A2 + Dall.Duu0{itype}'*(P.^2)*Dall.Duu0{itype}/Ntype_dxdy ...
            + Dall.Dvv0{itype}'*(P.^2)*Dall.Dvv0{itype}/Ntype_dxdy ...
            + 2*Dall.Duv0{itype}'*(P.^2)*Dall.Duv0{itype}/Ntype_dxdy;
    end
end

end

