function B = sfs_make_matrixB(omega2d)
% Make matrix B (boundary condition)
%
% -- Input
% omega2d : Object region. [Nx, Ny]
%
% -- Output
% B : Matrix for boundary condition. [Nomega, Nomega]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% preparation
omega = find(omega2d==1);
Nomega = numel(omega);
% boundary region
domega2d = 1 - double( conv2(1-omega2d,[0 1 0;1 1 1;0 1 0],'same') > 0);
domega2d = omega2d - domega2d;
domega = find(domega2d(omega)==1);
Ndomega = numel(domega);
% average coordinate values in boundary region
Nx = size(omega2d,1);
Ny = size(omega2d,2);
[meshY,meshX] = meshgrid(1:Ny,1:Nx);
xcm = sum(sum(domega2d.*meshX))/Ndomega;
ycm = sum(sum(domega2d.*meshY))/Ndomega;

%% create sparse matrix
b0i = ones(Ndomega^2,1);
b0j = ones(Ndomega^2,1);
b0s = zeros(Ndomega^2,1);
b1i = ones(Ndomega^2,1);
b1j = ones(Ndomega^2,1);
b1s = zeros(Ndomega^2,1);
bk = 0;
for i = 1 : Ndomega
    [ix,iy] = ind2sub([Nx,Ny],omega(domega(i)));
    for j = 1 : Ndomega
        [jx,jy] = ind2sub([Nx,Ny],omega(domega(j)));
        bk = bk + 1;
        b0i(bk) = domega(i);
        b0j(bk) = domega(j);
        b0s(bk) = 1/Ndomega^2;
        b1i(bk) = domega(i);
        b1j(bk) = domega(j);
        b1s(bk) = 1/Ndomega^2*(ix-xcm)*(jx-xcm) ...
            + 1/Ndomega^2*(iy-ycm)*(jy-ycm);
    end
end
B0 = sparse(b0i,b0j,b0s,Nomega,Nomega);
B1 = sparse(b1i,b1j,b1s,Nomega,Nomega);
B = B0 + B1;

end

