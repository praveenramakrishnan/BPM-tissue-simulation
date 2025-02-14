% Parameters 
lambda = 920e-9;

length_along_x = 10*lambda;
length_along_y = 10*lambda;
length_along_z = 5*lambda;

distance_to_focus = length_along_z;

delta.x = lambda/10;
delta.y = lambda/10;
delta.z = lambda/10;

I = round(length_along_x/delta.x);
J = round(length_along_y/delta.y);
K = round(length_along_z/delta.z);

K_focus = round(distance_to_focus/delta.z);

illorigin = [round(I/2), round(J/2), K_focus];
z_launch = 0;

% Number of slices into which region is divided. 
% This reduces memory requirement by loading refractive index data in smaller chunks.
num_region_slices = 2;

save(append(mfilename, '.mat'));
