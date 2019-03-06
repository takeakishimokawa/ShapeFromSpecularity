function parm_mfa1 = sfs_mfa1_beta(Ra, Duu, Dvv, H, parm_mfa1)
% Determine initial beta value for mean field algorithm
%
% -- Input
% Ra : Cholesky factorization of inverse of matrix A. Ra'*Ra=inv(A). [Nomega, Nomega]
% Duu : Matrix for differential oparator d^2/du^2. [Nomegain, Nomega]
% Dvv : Matrix for differential oparator d^2/dv^2. [Nomegain, Nomega]
% H : Diagonal matrix with elements (1-aniso). [Nomegain, Nomegain]
% parm_mfa1 : Parameters for mean field algorithm
%
% -- Output
% parm_mfa1 : Parameters for mean field algorithm
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

G = [Dvv' Duu'*H]; % J = G'*Ra'*Ra*G
parm_mfa1.beta = parm_mfa1.beta0/abs(eigs(Ra*(G*G')*Ra',1)); % eigs(J,1) = eigs((Ra*G)'*(Ra*G),1) = eigs((Ra*G)*(Ra*G)',1)

end

