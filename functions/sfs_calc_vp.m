function VP = sfs_calc_vp(im, omega2d, omega2d_orig, pyramidlevel, blurlevel)
% Calculate vertical polarity of intensity gradient
%
% -- Input
% im : raw image
% omega2d : Object region. [Nx, Ny]
% omega2d_orig : Object region of original image resolution
% pyramidlevel : Level of pyramid for filter response. Resolution is (original resolution)/2^(pyramidlevel).
% blurlevel : Filter response is downsampled to 1/2^blurlevel
%
% -- Output
% VP : Vertical polarity of intensity gradient. [Nx, Ny]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

%% object region of pyramidlevel
omega2d_pyr = double( blurDn(omega2d_orig,pyramidlevel) >= 0.25 );
% unreliable region near the boundary
domega2d_pyr = 1 - double( conv2(1-omega2d_pyr,ones(12),'same') > 0);
domega2d_pyr = omega2d_pyr - domega2d_pyr;

%% filter responce of vertical orientation
theta = pi/2; % vertical orientation
[pyr,pind,steermtx,harmonics] = buildSpyr(im,'auto','sp1Filters');
[lev,lind] = spyrLev(pyr,pind,pyramidlevel+1);
Nx_pyr = lind(1,1);
Ny_pyr = lind(1,2);
lev2 = reshape(lev,Nx_pyr*Ny_pyr,2);
P = reshape( steer(lev2, theta, harmonics, steermtx), Nx_pyr, Ny_pyr );
P(omega2d_pyr==0) = 0;
P(domega2d_pyr>0) = 0;
% downsample
Pdn = blurDn(P,blurlevel);
% Gaussian smoothing
mu = 9;
sigma = 4;
[meshY meshX] = meshgrid(1:(2*mu-1),1:(2*mu-1));
kernel = exp(-((meshX-mu).^2+(meshY-mu).^2)/(2*sigma^2));
kernel = kernel/sum(kernel(:));
Pdn = conv2(Pdn,kernel,'same');

%% obtain vertical polarity
VP = -sign(Pdn);
VP(omega2d==0) = 0;

end

