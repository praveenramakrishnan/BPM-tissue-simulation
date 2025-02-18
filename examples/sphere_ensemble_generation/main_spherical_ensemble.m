addpath('../../src/utils_simulation_setup/utils_sphere_code');
clear all; close all;
recompute_ensemble = true;
lambda = 920e-9;
length_along_x = 200*lambda;
length_along_y = 200*lambda;
length_along_z = 200*lambda;

delta = lambda/10;
x_grid_min = -length_along_x/2;
x_grid_max = length_along_x/2;
y_grid_min = -length_along_y/2;
y_grid_max = length_along_y/2;
z_grid_min = -length_along_z;
z_grid_max = 0;

radius_sphere_min = 0.5e-6;
radius_sphere_max = 0.5e-6;

% Concentration of cube is calculated using the web too at
%https://omlc.org/cgi-bin/mie_angles.cgi?diameter=10&lambda_vac=1.300&nr_sphere=1.39&ni_sphere=0&n_medium=1.3333&n_angles=100&density=0.00005
microns_cube = (1e-6)^3;
concentration_of_sphere = 9.4e-3/microns_cube;
volume_of_domain = length_along_x*length_along_y*length_along_z;
num_spheres = round(concentration_of_sphere*volume_of_domain);

filename_sphere_ensemble = 'sphere_ensemble.mat';
if recompute_ensemble
    [x_sphere_center, radius_sphere] = generate_ensemble( ...
        x_grid_min, x_grid_max, y_grid_min, y_grid_max, z_grid_min, z_grid_max, ...
        radius_sphere_min, radius_sphere_max, num_spheres);
    save(filename_sphere_ensemble, 'x_sphere_center', 'radius_sphere');
else
    load(filename_sphere_ensemble);
end

figure(1);clf;
plot3(x_sphere_center(1,:) ,x_sphere_center(2,:), x_sphere_center(3,:),'.');
axis equal;
