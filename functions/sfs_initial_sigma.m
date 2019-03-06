function [smax, smin] = sfs_initial_sigma(OF, aniso, VP, CCS, omegain)
% Calculate initial values of surface second derivative signs (smax and smin)
%
% -- Input
% OF : Image orientation. [Nx, Ny]
% aniso : Image anisotropy. [Nx, Ny]
% VP : Vertical polarity of intensity gradient. [Nx, Ny]
% CCS : Contour's curvature sign. [Nx, Ny]
% omegain : Index of omegain2d==1. [Nomegain, 1]
%
% -- Output
% smax : Initial value of large surface second derivative sign. [Nomegain, 1]
% smin : Initial value of small surface second derivative sign. [Nomegain, 1]
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

% divide object region into two regions
omegain_a = ( sin(OF(omegain)).^2 >= (1-aniso(omegain)).*cos(OF(omegain)).^2 );
omegain_b = ( sin(OF(omegain)).^2 < (1-aniso(omegain)).*cos(OF(omegain)).^2 );

% divided vertical polarity
smax = VP(omegain);
smax(omegain_b) = 0;
smin = VP(omegain);
smin(omegain_a) = 0;

% initial value near boundary is determined by contour's curvature sign
omegain_ccs = find(abs(CCS(omegain))==1);
smax(omegain_ccs) = 1;
smin(omegain_ccs) = CCS(omegain(omegain_ccs));

end

