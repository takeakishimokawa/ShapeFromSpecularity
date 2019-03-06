function CCS = sfs_calc_ccs(omega2d, omega2d_orig, downlevel)
% Calculate contour's curvature sign
%
% -- Input
% omega2d : Object region. [Nx, Ny]
% omega2d_orig : Object region of original image resolution
% downlevel : Level of downsampling. 
%
% -- Output
% CCS : Contour's curvature sign. [Nx, Ny]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% preparation
% convolution kernel
Nx_orig = size(omega2d_orig,1);
Ny_orig = size(omega2d_orig,2);
Nim_orig = min([Nx_orig,Ny_orig]);
kernel_a = [0 1 0;1 1 1;0 1 0]; % adjacency kernel
radius_c = Nim_orig/8;
radius_s = Nim_orig/64;
[meshY meshX] = meshgrid(1:(2*radius_c+1),1:(2*radius_c+1));
kernel_c = double( (meshX-radius_c-1).^2+(meshY-radius_c-1).^2 <= radius_c^2 ); % circular kernel to measure curvature sign
[meshY meshX] = meshgrid(1:(2*radius_s+1),1:(2*radius_s+1));
kernel_s = double( (meshX-radius_s-1).^2+(meshY-radius_s-1).^2 <= radius_s^2 ); % circular kernel for smoothing
% boundary points
inner_boundary =  1 - double( conv2(1-omega2d_orig,kernel_a,'same') > 0);
inner_boundary = omega2d_orig - inner_boundary;
outer_boundary = double( conv2(omega2d_orig,kernel_a,'same') > 0);
outer_boundary = outer_boundary - omega2d_orig;
% extend object region to avoid errors that occur when boundary points are near image edge
omega2d_orig_extend = zeros(Nx_orig+2*radius_c,Ny_orig+2*radius_c);
omega2d_orig_extend(radius_c+(1:Nx_orig),radius_c+(1:Ny_orig)) = omega2d_orig;

%% calculate contour's curvature sign
CCS = conv2(1-2*omega2d_orig_extend,kernel_c,'same');
CCS = CCS(radius_c+(1:Nx_orig),radius_c+(1:Ny_orig));
CCS(inner_boundary==0&outer_boundary==0) = 0;
CCS = sign(CCS);
% smoothing
CCS = sign(conv2(CCS,kernel_s,'same'));
CCS(omega2d_orig==0) = 0;
% downsampling
CCSregion = blurDn(abs(CCS),downlevel);
CCSregion(CCSregion<0.25) = 0;
CCS = blurDn(CCS,downlevel);
CCS = sign(CCS);
CCS(omega2d==0) = 0;
CCS(CCSregion==0) = 0;

end

