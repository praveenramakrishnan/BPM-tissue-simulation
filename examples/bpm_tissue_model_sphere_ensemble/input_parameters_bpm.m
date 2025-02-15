% Parameters 
lambda = 920e-9; % wavelength

% Dimensions of the domain
length_along_x = 50*lambda;
length_along_y = 50*lambda;
length_along_z = 5*lambda;

% Distance from start of the domain to the focus
distance_to_focus = length_along_z;

% Grid discretisation size
delta.x = lambda/10;
delta.y = lambda/10;
delta.z = lambda/10;

% Number of points along each axis
I = round(length_along_x/delta.x);
J = round(length_along_y/delta.y);
K = round(length_along_z/delta.z);

% The number of cells from the start to focus
K_focus = round(distance_to_focus/delta.z);

% Local of the origin of the domain in terms of number of cells
illorigin = [round(I/2), round(J/2), K_focus];

% Additional shift in z axis
z_launch = 0;

% Number of planes used in BPM computation
num_bpm_planes = 150;

% Number of slices into which region is divided.
% This reduces memory requirement by loading refractive index data in smaller chunks.
num_region_slices = 2;

save(append(mfilename, '.mat'));
