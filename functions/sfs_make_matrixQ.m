function Q = sfs_make_matrixQ(Dall, Nomega, Nomegain)
% Make matrix Q used in gradient descent
%
% -- Input
% Dall : Full data for making matrix D
% Nomega : number of pixels in omega2d
% Nomegain : number of pixels in omegain2d
%
% -- Output
% Q : Matries used in gradient descent. {Nomegain, 1}
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

Q = cell(Nomegain,1);
Ntype_dxdy = size(Dall.Duu0,1);
for i = 1 : Nomegain
    tmp = sparse(1,1,0,Nomega,Nomega);
    Ii = sparse(i,i,1,Nomegain,Nomegain);
    for itype = 1 : Ntype_dxdy
        tmp = tmp + Dall.Duu0{itype}'*Ii*Dall.Duu0{itype}/Ntype_dxdy ...
            + Dall.Dvv0{itype}'*Ii*Dall.Dvv0{itype}/Ntype_dxdy ...
            + 2*Dall.Duv0{itype}'*Ii*Dall.Duv0{itype}/Ntype_dxdy;
    end
    Q{i} = tmp;
end

end

