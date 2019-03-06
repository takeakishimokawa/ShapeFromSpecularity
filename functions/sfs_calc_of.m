function [OF, aniso, OFvis] = sfs_calc_of(im, omega2d, Ndir, pyramidlevel, blurlevel)
% Calculate orientation field
%
% -- Input
% im : raw image
% omega2d : Object region. [Nx, Ny]
% Ndir : Angle resolution of image orientation (0-180deg)
% pyramidlevel : Level of pyramid we use for filter response. Resolution is (original resolution)/2^(pyramidlevel)
% blurlevel : Magnitude of filter response is downsampled to 1/2^blurlevel
%
% -- Output
% OF : Image orientation. [Nx, Ny]
% aniso : Image anisotropy. [Nx, Ny]
% OFvis : Orientation field data for visualization. Hue represents image
% orientation. Saturation represents image anisotropy. [Nx, Ny]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% filter responce in each orientation
theta = pi * [0:Ndir-1]'/Ndir;
[pyr,indices,steermtx,harmonics] = buildSpyr(im,'auto','sp1Filters');
[lev,lind] = spyrLev(pyr,indices,pyramidlevel+1);
Nx_pyr = lind(1,1);
Ny_pyr = lind(1,2);
lev2 = reshape(lev,Nx_pyr*Ny_pyr,2);
P = zeros(Ndir,Nx_pyr,Ny_pyr);
for i = 1 : Ndir
    P(i,:,:) = reshape( steer(lev2,theta(i),harmonics,steermtx), Nx_pyr, Ny_pyr );
end

%% magnitude of filter response
P2 = P.^2;
% downsample
Nx = Nx_pyr/2^blurlevel;
Ny = Ny_pyr/2^blurlevel;
P2dn = zeros(Ndir,Nx,Ny);
for i = 1 : Ndir
    P2dn(i,:,:) = blurDn(squeeze(P2(i,:,:)),blurlevel);
end
% smoothing
for i = 1 : Ndir
    P2dn(i,:,:) = conv2(squeeze(P2dn(i,:,:)),ones(3),'same');
end

%% obtain image orientation and image anisotropy
OF = zeros(Nx,Ny);
aniso = zeros(Nx,Ny);
OFvis = uint8(256*ones(Nx,Ny,3)); % for visualization
color_rgb = 1-hsv(Ndir);
for i = 1 : Nx
    for j = 1 : Ny
        Pij = P2dn(:,i,j);
        [Cmax Imax] = max(Pij);
        [Cmin Imin] = min(Pij);
        if omega2d(i,j)>0
            aniso(i,j) = 1-sqrt(Cmin/Cmax);
            g = mod(Imax-1+0,Ndir)+1;
            OF(i,j) = theta(g);
            OFvis(i,j,1) = uint8(256*(1-color_rgb(g,1)*aniso(i,j)));
            OFvis(i,j,2) = uint8(256*(1-color_rgb(g,2)*aniso(i,j)));
            OFvis(i,j,3) = uint8(256*(1-color_rgb(g,3)*aniso(i,j)));
        end
    end
end

end

