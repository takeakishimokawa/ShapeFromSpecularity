function [r_g, r_li, r_g0, r_li0] = sfs_depth_corr(z2d, true_z2d, omega2d)
% Calculate depth correlation between recovered and true shape
%
% -- Input
% z2d : Recovered 3D shape. [Nx, Ny]
% true_z2d : True 3D shape. [Nx, Ny]
% omega2d : Object region. [Nx, Ny]
%
% -- Output
% r_g : Global depth correlation
% r_li : Local interior depth correlation
% r_g0 : Global depth correlation without slant calibration
% r_li0 : Local interior depth correlation without slant calibration
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% preparation
Nx = size(omega2d,1);
Ny = size(omega2d,2);
Nim = min([Nx,Ny]);
omega = find(omega2d==1);
z = z2d(omega);
true_z = true_z2d(omega);
% boundary region
domega2d = 1 - double( conv2(1-omega2d,[0 1 0;1 1 1;0 1 0],'same') > 0);
domega2d = omega2d - domega2d;
domega = find(domega2d(omega)==1);
Ndomega = numel(domega);
% slant calibration
[meshY,meshX] = meshgrid(1:Ny,1:Nx);
xcm = sum(sum(domega2d.*meshX))/Ndomega;
ycm = sum(sum(domega2d.*meshY))/Ndomega;
A_slant = [sum(sum(domega2d.*(meshX-xcm).^2)), sum(sum(domega2d.*(meshX-xcm).*(meshY-ycm))); ...
    sum(sum(domega2d.*(meshX-xcm).*(meshY-ycm))), sum(sum(domega2d.*(meshY-ycm).^2))];
b_slant = - [sum(sum(domega2d.*(meshX-xcm).*true_z2d)); sum(sum(domega2d.*(meshY-ycm).*true_z2d))];
x_slant = A_slant\b_slant;
true_z2d_noslant = true_z2d + x_slant(1)*(meshX-xcm) + x_slant(2)*(meshY-ycm);
true_z_noslant = true_z2d_noslant(omega);

%% global depth correlation
r_g = corr(z,true_z_noslant);
r_g0 = corr(z,true_z);

%% local interior correlation
% removal of the area near the boundary
radius_rb = 24*Nim/256;
[meshY meshX] = meshgrid(1:(2*radius_rb+1),1:(2*radius_rb+1));
kernel_rb = double( (meshX-radius_rb-1).^2+(meshY-radius_rb-1).^2 <= radius_rb^2 );
omega2d_rb = double( conv2(1-omega2d,kernel_rb,'same') == 0 );
% circular kernel for interior correlation
Ngrid = 8;
radius_g = Nim/Ngrid;
[meshY meshX] = meshgrid(1:(2*radius_g+1),1:(2*radius_g+1));
kernel_g = double( (meshX-radius_g-1).^2+(meshY-radius_g-1).^2 <= radius_g^2 );
% local interior correlation of each grid
R_li_noslant = NaN(Ngrid-1,Ngrid-1);
R_li = NaN(Ngrid-1,Ngrid-1);
for ix = 1 : (Ngrid-1)
    for iy = 1 : (Ngrid-1)
        local_center = zeros(Nx,Ny);
        local_center(radius_g*ix,radius_g*iy) = 1;
        local_region2d = conv2(local_center,kernel_g,'same');
        overlap_local_region2d = omega2d_rb.*local_region2d;
        overlap_local_region = find(overlap_local_region2d==1);
        if sum(overlap_local_region2d(:)) > 0.5*sum(local_region2d(:))
            R_li_noslant(ix,iy) = corr(z2d(overlap_local_region),true_z2d_noslant(overlap_local_region));
            R_li(ix,iy) = corr(z2d(overlap_local_region),true_z2d(overlap_local_region));
        end
    end
end
r_li = nanmean(R_li_noslant(:));
r_li0 = nanmean(R_li(:));

end
