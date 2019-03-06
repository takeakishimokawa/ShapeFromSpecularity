function [smax, smin] = sfs_mfa2(smax, smin, ka, Ra2, Duu, Dvv, H, CCS, omegain, parm_mfa2)
% Mean field algorithm for second cost function
%
% -- Input
% smax : Large surface second derivative sign. [Nomegain, 1]
% smin : Small surface second derivative sign. [Nomegain, 1]
% ka : Curvature magnitudes. [Nomegain, 1]
% Ra2 : Cholesky factorization of inverse of matrix A2. Ra2'*Ra2=inv(A2). [Nomega, Nomega]
% Duu : Matrix for differential oparator d^2/du^2. [Nomegain, Nomega]
% Dvv : Matrix for differential oparator d^2/dv^2. [Nomegain, Nomega]
% H : Diagonal matrix with elements (1-aniso). [Nomegain, Nomegain]
% CCS : Contour's curvature sign. [Nx, Ny]
% omegain : Index of omegain2d==1. [Nomegain, 1]
% parm_mfa2 : Parameters for mean field algorithm
%
% -- Output
% smax : Optimized large principal curvature sign. [Nomegain, 1]
% smin : Optimized small principal curvature sign. [Nomegain, 1]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

sall = [smax; smin];
Nomegain = numel(omegain);
P = spdiags(1./ka,0,Nomegain,Nomegain);
G2 = [Dvv'*P Duu'*P*H]; % J2 = G2'*Ra2'*Ra2*G2
hall = [abs(CCS(omegain)); CCS(omegain)];

beta2 = parm_mfa2.beta0/abs(eigs(Ra2*(G2*G2')*Ra2',1)); % eigs(J2,1) = eigs((Ra2*G2)'*(Ra2*G2),1) = eigs((Ra2*G2)*(Ra2*G2)',1)
for repeat = 1 : parm_mfa2.N
    beta2 = beta2 * parm_mfa2.beta_gain;
    sall = tanh(beta2*(G2'*(Ra2'*(Ra2*(G2*sall)))+hall));
end
sall = sign(sall);
smax = sall(1:Nomegain);
smin = sall((Nomegain+1):end);

end

