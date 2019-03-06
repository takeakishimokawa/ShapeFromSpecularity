function [Duu, Duv, Dvv, Duu2, Duv2, Dvv2, Dall] = sfs_make_matrixD(OF, omega2d)
% Make matrix D (second order differential oparator)
%
% -- Input
% OF : Image orientation. [Nx, Ny]
% omega2d : Object region. [Nx, Ny]
%
% -- Output
% Duu : Matrix for differential oparator d^2/du^2. [Nomegain, Nomega]
% Duv : Matrix for differential oparator d^2/dudv. [Nomegain, Nomega]
% Dvv : Matrix for differential oparator d^2/dv^2. [Nomegain, Nomega]
% Duu2 : Duu'*Duu. [Nomega, Nomega]
% Duv2 : Duv'*Duv. [Nomega, Nomega]
% Dvv2 : Dvv'*Dvv. [Nomega, Nomega]
% Dall : Full data for making matrix D
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% preparation
Nx = size(omega2d,1);
Ny = size(omega2d,2);
omega = find(omega2d==1);
Nomega = numel(omega);
omegain2d = double( conv2(omega2d,ones(3),'same') >= 9 ); % omegain2d is the interior of omega2d (surrounded by 9 pixels). Second derivative values exist here.
omegain = find(omegain2d==1);
Nomegain = numel(omegain);
% 9 surrounding pixels are used for second order differential operation.
xy_surr = [ ... % representing coordinate value of surrounding pixel
     0  0; % (i  ,j  ) itself
     1  0; % (i+1,j  ) right
     1  1; % (i+1,j+1) upper right
     0  1; % (i  ,j+1) upper
    -1  1; % (i-1,j+1) upper left
    -1  0; % (i-1,j  ) left
    -1 -1; % (i-1,j-1) lower left
     0 -1; % (i  ,j-1) lower
     1 -1];% (i+1,j-1) lower right
Nsurr = size(xy_surr,1); % = 9
i_omega = cell(Nsurr,1);
% d^2/dxdy is calculated in 4 ways and averaged so that matrix A does not become singular.
type_dxdy = [ ...
    1 -1 1 -1 0 0 0 0 0; % f(i+1,j+1) - f(i+1,j) - f(i,j+1) + f(i,j)
    -1 0 0 1 -1 1 0 0 0; % f(i,j+1) - f(i,j) - f(i-1,j+1) + f(i-1,j)
    1 0 0 0 0 -1 1 -1 0; % f(i,j) - f(i,j-1) - f(i-1,j) + f(i-1,j-1)
    -1 1 0 0 0 0 0 1 -1]'; % f(i+1,j) - f(i+1,j-1) - f(i,j) + f(i,j-1)
Ntype_dxdy = size(type_dxdy,2); % = 4
% d^2/dx^2 is calculated by [-2 1 0 0 0 1 0 0 0]'
% d^2/dy^2 is calculated by [-2 0 0 1 0 0 0 1 0]'

%% create sparse matrix
% for Duu
duui = ones(Nomegain*Nsurr,Ntype_dxdy);
duuj = ones(Nomegain*Nsurr,Ntype_dxdy);
duus = zeros(Nomegain*Nsurr,Ntype_dxdy);
% for Duv
duvi = ones(Nomegain*Nsurr,Ntype_dxdy);
duvj = ones(Nomegain*Nsurr,Ntype_dxdy);
duvs = zeros(Nomegain*Nsurr,Ntype_dxdy);
% for Dvv
dvvi = ones(Nomegain*Nsurr,Ntype_dxdy);
dvvj = ones(Nomegain*Nsurr,Ntype_dxdy);
dvvs = zeros(Nomegain*Nsurr,Ntype_dxdy);
dk = 0;
for i_omegain = 1 : Nomegain
    [ix0,iy0] = ind2sub([Nx,Ny],omegain(i_omegain));
    for isurr = 1 : Nsurr
        ix = ix0 + xy_surr(isurr,1);
        iy = iy0 + xy_surr(isurr,2);
        i_omega{isurr} = find(omega==ix+(iy-1)*Nx);
    end
    theta = OF(ix0,iy0);
    ct = cos(theta);
    st = sin(theta);
    % d^2/du^2 = d^2/dx^2 ct^2 + 2*d^2/dxdy ct*st + d^2/dy^2 st^2
    duu0 = [ ...  % duu0 = d^2/dx^2 ct^2 + d^2/dy^2 st^2
        - 2*ct^2 - 2*st^2;
        ct^2;
        0;
        st^2;
        0;
        ct^2;
        0;
        st^2;
        0];
    % d^2/dudv = -d^2/dx^2 ct*st + d^2/dxdy (ct^2-st^2) + d^2/dy^2 ct*st
    duv0 = [ ... % duv0 = -d^2/dx^2 ct*st + d^2/dy^2 ct*st
        2*ct*st - 2*ct*st;
        -ct*st;
        0;
        ct*st;
        0;
        -ct*st;
        0;
        ct*st;
        0];
    % d^2/dv^2 = d^2/dx^2 st^2 - 2*d^2/dxdy ct*st + d^2/dy^2 ct^2
    dvv0 = [ ... % dvv0 = d^2/dx^2 st^2 + d^2/dy^2 ct^2
        - 2*st^2 - 2*ct^2;
        st^2;
        0;
        ct^2;
        0;
        st^2;
        0;
        ct^2;
        0];
    for isurr = 1 : Nsurr
        dk = dk + 1;
        for itype = 1 : Ntype_dxdy
            duu = duu0 + type_dxdy(:,itype)*2*ct*st; % add 2*d^2/dxdy ct*st
            duv = duv0 + type_dxdy(:,itype)*(ct^2-st^2); % add d^2/dxdy (ct^2-st^2)
            dvv = dvv0 - type_dxdy(:,itype)*2*ct*st; % add -2*d^2/dxdy ct*st
            duui(dk,itype) = i_omegain;
            duuj(dk,itype) = i_omega{isurr};
            duus(dk,itype) = duu(isurr);
            duvi(dk,itype) = i_omegain;
            duvj(dk,itype) = i_omega{isurr};
            duvs(dk,itype) = duv(isurr);
            dvvi(dk,itype) = i_omegain;
            dvvj(dk,itype) = i_omega{isurr};
            dvvs(dk,itype) = dvv(isurr);
        end
    end
end
Duu0 = cell(Ntype_dxdy,1);
Duv0 = cell(Ntype_dxdy,1);
Dvv0 = cell(Ntype_dxdy,1);
for itype = 1 : Ntype_dxdy
    Duu0{itype} = sparse(duui(:,itype),duuj(:,itype),duus(:,itype),Nomegain,Nomega);
    Duv0{itype} = sparse(duvi(:,itype),duvj(:,itype),duvs(:,itype),Nomegain,Nomega);
    Dvv0{itype} = sparse(dvvi(:,itype),dvvj(:,itype),dvvs(:,itype),Nomegain,Nomega);
end

%% averaging for stabilization
Duu = sparse(Nomegain,Nomega);
Duv = sparse(Nomegain,Nomega);
Dvv = sparse(Nomegain,Nomega);
Duu2 = sparse(Nomega,Nomega);
Duv2 = sparse(Nomega,Nomega);
Dvv2 = sparse(Nomega,Nomega);
for itype = 1 : Ntype_dxdy
    Duu = Duu + Duu0{itype}/Ntype_dxdy;
    Duv = Duv + Duv0{itype}/Ntype_dxdy;
    Dvv = Dvv + Dvv0{itype}/Ntype_dxdy;
    Duu2 = Duu2 + Duu0{itype}'*Duu0{itype}/Ntype_dxdy;
    Duv2 = Duv2 + Duv0{itype}'*Duv0{itype}/Ntype_dxdy;
    Dvv2 = Dvv2 + Dvv0{itype}'*Dvv0{itype}/Ntype_dxdy;
end

%%
Dall.Duu0 = Duu0;
Dall.Duv0 = Duv0;
Dall.Dvv0 = Dvv0;

end

