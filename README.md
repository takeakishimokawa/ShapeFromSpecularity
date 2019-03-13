# ShapeFromSpecularity
MATLAB code for estimating 3D shape from a single specular image proposed in *Shimokawa et al. 2019*.
* [Takeaki Shimokawa, Akiko Nishio, Masa-aki Sato, Mitsuo Kawato, and Hidehiko Komatsu, "Computational model for human 3D shape perception from a single specular image," Frontiers in Computational Neuroscience 13 (2019) 10](https://doi.org/10.3389/fncom.2019.00010).

Contact [Takeaki Shimokawa](https://sites.google.com/view/takeakishimokawa/) if you have any questions.

## Setting up and running a demo code
1. Clone or download our MATLAB code.
2. Our MATLAB code requires an external MATLAB toolbox. Download matlabPyrTools from
https://github.com/LabForComputationalVision/matlabPyrTools
to *./external directory*.
3. Run *sfs_demo.m*

#### Expected run time
If the resultant resolution of the 3D shape recovery is set to 128 x 128 (downlevel=3), the run time takes about ten minutes (desktop CPU, 3 GHz) and requires 4+GB RAM. If the resultant resolution is set to 256 x 256 (downlevel=2), it takes about three hours and requires 16+GB RAM.

#### Troubleshooting
If you have MEX-file error, recompile the MEX-files of the matlabPyrTools as follows:
1. Move to *./external/matlabPyrTools-master/MEX/*
2. Run *compilePyrTools.m*
3. Copy the recompiled MEX-files in *./external/matlabPyrTools-master/MEX/* to *./external/matlabPyrTools-master/*

## Description of functions and images
* *sfs_demo.m* : Demo function for estimating 3D shape from an image.
* *functions/sfs_calc_of.m* : This function calculates the orientation field from an image.
* *functions/sfs_calc_vp.m* : This function calculates the vertical polarity of the intensity gradient from an image.
* *functions/sfs_calc_ccs.m* : This function calculates the contour's curvature sign from an object region image.
* *functions/sfs_initial_sigma.m* : This function calculates the initial values of the surface second derivative signs.
* *functions/sfs_mfa1.m* : This function implements the mean field algorithm based on the first cost function.
* *functions/sfs_mfa1_beta.m* : This function determines the initial beta value required for *sfs_mfa1.m*.
* *functions/sfs_mfa2.m* : This function implements the mean field algorithm based on the second cost function.
* *functions/sfs_ocm.m* : This function optimizes the second derivative magnitudes and signs based on the first cost function.
* *functions/sfs_gd.m* : This function implements a gradient descent based on the second cost function.
* *functions/sfs_make_matrixB.m* : This function makes matrix B (boundary condition).
* *functions/sfs_make_matrixD.m* : This function makes matrix D (second order differential operator).
* *functions/sfs_make_matrixQ.m* : This function makes matrix Q, which is required for *sfs_gd.m*.
* *functions/sfs_depth_corr.m* : This function calculates the depth correlation between the recovered and ground-truth shapes.
* *data/image/image1-26.tif* : Glossy surface images (1-12), mirrored surface images (13-24), and images used for psychophysical experiment in *Shimokawa et al. 2019* (25-26). All of them had 1024 x 1024 pixel resolution.
* *data/image/objectregion/image1-26ob.tif* : Images indicating the object region (silhouette).
* *data/shape/depth1-26.mat* : Ground-truth 3D shapes for evaluating the shape recovery performance (256 x 256 or 128 x 128 resolution).
