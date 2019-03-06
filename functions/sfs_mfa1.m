function [smax, smin] = sfs_mfa1(smax, smin, Ra, Duu, Dvv, H, CCS, omegain, parm_mfa1)
% Mean field algorithm for optimizing surface second derivative signs (smax and smin)
%
% -- Input
% smax : Large surface second derivative sign. [Nomegain, 1]
% smin : Small surface second derivative sign. [Nomegain, 1]
% Ra : Cholesky factorization of inverse of matrix A. Ra'*Ra=inv(A). [Nomega, Nomega]
% Duu : Matrix for differential oparator d^2/du^2. [Nomegain, Nomega]
% Dvv : Matrix for differential oparator d^2/dv^2. [Nomegain, Nomega]
% H : Diagonal matrix with elements (1-aniso). [Nomegain, Nomegain]
% CCS : Contour's curvature sign. [Nx, Ny]
% omegain : Index of omegain2d==1. [Nomegain, 1]
% parm_mfa1 : Parameters for mean field algorithm
%
% -- Output
% smax : Optimized large principal curvature sign. [Nomegain, 1]
% smin : Optimized small principal curvature sign. [Nomegain, 1]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

sall = [smax; smin];
G = [Dvv' Duu'*H]; % J = G'*Ra'*Ra*G
hall = [abs(CCS(omegain)); CCS(omegain)];
Nomegain = numel(omegain);

beta = parm_mfa1.beta;
for repeat = 1 : parm_mfa1.N
    beta = beta * parm_mfa1.beta_gain;
    sall = tanh(beta*(G'*(Ra'*(Ra*(G*sall)))+hall));
end
sall = sign(sall);
smax = sall(1:Nomegain);
smin = sall((Nomegain+1):end);

end

