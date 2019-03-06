% Demo program
% MATLAB code for estimating 3D shape from a single specular image.
%
% Copyright (C) 2019, Takeaki Shimokawa, ATR.

close all
clear
tic

%% 0. Add path
addpath(genpath('./external'))
addpath(genpath('./functions'))

%% 1. Load an specular image
imageID = 2; % glossy surface: 1-12, mirrored surface: 13-24, images used for psychophysical experiment: 25-26

im = imread(sprintf('./data/image/image%d.tif',imageID));
figure('Name','Specular image'),image(im),axis square,axis off % raw image
im = mean(im,3);

omega2d_orig = imread(sprintf('./data/image/objectregion/image%dor.tif',imageID)); % omega: object region
omega2d_orig = 1-sign(mean(omega2d_orig,3));

%% 2. Calculate orientation field, vertical polarity, and contour's curvature sign 
downlevel = 3;%2; % Level of downsampling. Image size will be 1/2^downlevel.

% downsample object region
omega2d = double( blurDn(omega2d_orig,downlevel) >= 0.25 );
% opening to avoid problems in differential operation
omega2d = double( conv2(omega2d,ones(3),'same') >= 9 ); % erosion
omega2d = double( conv2(omega2d,ones(3),'same') > 0 ); % dilation
% indices of object region etc.
omega = find(omega2d==1);
Nomega = numel(omega);
omegain2d = double( conv2(omega2d,ones(3),'same') >= 9 ); % omegain2d is the interior of omega2d. Second derivative values exist here.
omegain = find(omegain2d==1);
Nomegain = numel(omegain);
Nx = size(omega2d,1);
Ny = size(omega2d,2);

% orientation field
pyramidlevel = 0; % choose finest filter response
blurlevel = downlevel - pyramidlevel;
Ndir = 120; % angle resolution of image orientation (0-180deg)
[OF,aniso,OFvis] = sfs_calc_of(im,omega2d,Ndir,pyramidlevel,blurlevel);
figure('Name','Orientation field'),image(OFvis),axis square,axis off % visualize orientation field

% vertical polarity
pyramidlevel = min(downlevel,2); % choose relatively low-frequency filter response
blurlevel = downlevel - pyramidlevel;
VP = sfs_calc_vp(im,omega2d,omega2d_orig,pyramidlevel,blurlevel);
figure('Name','Vertical polarity of intensity gradient'),imagesc(VP),colormap(gray),axis square,axis off % visualize vertical polarity

% contour's curvature sign
CCS = sfs_calc_ccs(omega2d,omega2d_orig,downlevel);
% figure('Name','Contour''s curvature sign'),imagesc(CCS),colormap(gray),axis square,axis off % visualize contour's curvature sign

%% 3. Formulate cost function
B = sfs_make_matrixB(omega2d);
[Duu,Duv,Dvv,Duu2,Duv2,Dvv2,Dall] = sfs_make_matrixD(OF,omega2d);
A = Duu2 + Dvv2 + 2*Duv2 + B;
invA = inv(full(A));
Ra = chol(invA); % Ra'*Ra=inv(A)
H = spdiags(1-aniso(omegain),0,Nomegain,Nomegain);

%% 4. Minimize first cost function
Nloop = 10; % number of loop (= mean field algorithm & optimizaion of second derivative magnitudes)
% parameters for mean field algorithm
parm_mfa1.N = 100; % number of iteration
parm_mfa1.beta0 = 10;
parm_mfa1.beta_gain = 1.1;
parm_mfa1 = sfs_mfa1_beta(Ra,Duu,Dvv,H,parm_mfa1);
% parameters for optimizaion of second derivative magnitudes
parm_ocm.N = 10; % number of iteration

% initial sigma values
[smax,smin] = sfs_initial_sigma(OF,aniso,VP,CCS,omegain);
% smax2d=zeros(Nx,Ny);smax2d(omegain)=smax;figure('Name','Initial smax'),imagesc(smax2d),colormap(gray),axis square,axis off % visualize initial smax
% smin2d=zeros(Nx,Ny);smin2d(omegain)=smin;figure('Name','Initial smin'),imagesc(smin2d),colormap(gray),axis square,axis off % visualize initial smin

% loop of mean field algorithm & optimizaion of second derivative magnitudes
for loop = 1 : Nloop
    % mean field algorithm
    [smax,smin] = sfs_mfa1(smax,smin,Ra,Duu,Dvv,H,CCS,omegain,parm_mfa1);
%     smax2d=zeros(Nx,Ny);smax2d(omegain)=smax;figure('Name','Optimized smax'),imagesc(smax2d),colormap(gray),axis square,axis off % visualize optimized smax
%     smin2d=zeros(Nx,Ny);smin2d(omegain)=smin;figure('Name','Optimized smin'),imagesc(smin2d),colormap(gray),axis square,axis off % visualize optimized smin
    % optimizaion of second derivative magnitudes
    for repeat = 1 : parm_ocm.N
        [ka,smax,smin] = sfs_ocm(smax,smin,invA,Duu,Dvv,H,Nomegain);
    end
end

% intermediate solution
Ka = spdiags(ka,0,Nomegain,Nomegain);
z = - A \ ( Dvv'*Ka*smax + Duu'*Ka*H*smin );
z2d = max(z(:))*ones(Nx,Ny);
z2d(omega) = z;
% figure('Name','Intermediate solution'),imagesc(z2d),colormap(gray),hold on,contour(z2d,15,'k'),axis square,axis off % visualize intermediate solution
clearvars invA Ra

%% 5. Minimize second cost function
Nloop2 = 10; % number of loop (= mean field algorithm & gradient descent)
% parameters for gradient descent
parm_gd.N = 100; % number of iteration
parm_gd.gamma = 0.1; % step size
% parameters for mean field algorithm
parm_mfa2.N = 100; % number of iteration
parm_mfa2.beta0 = 10;
parm_mfa2.beta_gain = 1.1;

% initial values etc.
ka = ones(Nomegain,1);
A2 = A;
Q = sfs_make_matrixQ(Dall,Nomega,Nomegain);

% loop of gradient descent & mean field algorithm
for loop = 1 : Nloop2
    % gradient descent
    [ka,A2] = sfs_gd(ka,smax,smin,A2,B,Duu,Dvv,Dall,H,Q,Nomegain,parm_gd);
    if loop == Nloop2
        break
    end
    Ra2 = chol(inv(full(A2))); % Ra2'*Ra2=inv(A2)
    % mean field algorithm
    [smax,smin] = sfs_mfa2(smax,smin,ka,Ra2,Duu,Dvv,H,CCS,omegain,parm_mfa2);
%     smax2d=zeros(Nx,Ny);smax2d(omegain)=smax;figure('Name','Optimized smax'),imagesc(smax2d),colormap(gray),axis square,axis off % visualize optimized smax
%     smin2d=zeros(Nx,Ny);smin2d(omegain)=smin;figure('Name','Optimized smin'),imagesc(smin2d),colormap(gray),axis square,axis off % visualize optimized smin
end

% recovered 3D shape
P = spdiags(1./ka,0,Nomegain,Nomegain);
z = - A2 \ ( Dvv'*P*smax + Duu'*P*H*smin );
z2d = max(z(:))*ones(Nx,Ny);
z2d(omega) = z;
figure('Name','Recovered 3D shape'),imagesc(z2d),colormap(gray),hold on,contour(z2d,15,'k'),axis square,axis off % visualize recovered 3D shape

%% 6. Evaluate shape recovery performance
% for only downlevel=2 (256x256 resolution) or downlevel=3 (128x128 resolution)
load(sprintf('./data/shape/depth%d.mat',imageID))
if downlevel == 2
    [r_g,r_li] = sfs_depth_corr(z2d,true_z2d_256,omega2d_256)
%     true_z2d_256(omega2d_256<1) = max(true_z2d_256(:));
%     figure('Name','True 3D shape'),imagesc(true_z2d_256),colormap(gray),hold on,contour(true_z2d_256,15,'k'),axis square,axis off % visualize true 3D shape
elseif downlevel == 3
    [r_g,r_li] = sfs_depth_corr(z2d,true_z2d_128,omega2d_128)
%     true_z2d_128(omega2d_128<1) = max(true_z2d_128(:));
%     figure('Name','True 3D shape'),imagesc(true_z2d_128),colormap(gray),hold on,contour(true_z2d_128,15,'k'),axis square,axis off % visualize true 3D shape
end

toc
