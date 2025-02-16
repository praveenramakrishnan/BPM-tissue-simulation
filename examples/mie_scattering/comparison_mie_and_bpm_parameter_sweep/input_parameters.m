% Parameters 
lambda = 1300e-9; % wavelength

% Dimensions of the domain
length_along_x = 10*lambda;
length_along_y = 10*lambda;
length_along_z = 10*lambda;

% Grid discretisation size
delta.x = lambda/10;
delta.y = lambda/10;
delta.z = lambda/10;

% Number of points along each axis
I = round(length_along_x/delta.x);
J = round(length_along_y/delta.y);
K = round(length_along_z/delta.z);

% Local of the origin of the domain in terms of number of cells
illorigin = [round(I/2), round(J/2), round(K/2)];

% Additional shift in z axis
z_launch = 0;

% Number of planes used in BPM computation
num_bpm_planes = 150;

% Number of slices into which region is divided for BPM computations.
% This reduces memory requirement by loading refractive index data in smaller chunks.
num_region_slices = 1;

% Refractive index and sphere parameters
refractive_index_sphere_sweep_list = [1.375, 1.4];
refractive_index_background = 1.3333;
radius_sphere_sweep_list = [lambda, 2*lambda];

% Minimum number of terms in mie series calculation
num_terms_mie = 50;

save(append(mfilename, '.mat'));
