clear all; close all;
% Add path
addpath('../../../src/utils_mie_scattering/');

%Parameters
num_terms_mie = 50;
refractive_index_sphere = 1.3;
refractive_index_background = 1.0;
efield_amplitude_incident = 1.0;
radius_sphere = 10e-6;
wavelength = 1300e-9;

num_grid_points_x = 201;
num_grid_points_z = 201;
width_along_x = 30e-6;
length_along_z = 25e-6;

x_grid = linspace(-width_along_x/2, width_along_x/2, num_grid_points_x);
y_grid = 0;
z_grid = linspace(0, length_along_z, num_grid_points_z);

[x_grid_nd, y_grid_nd, z_grid_nd] = ndgrid(x_grid, y_grid, z_grid);
vertices_grid = [x_grid_nd(:), y_grid_nd(:), z_grid_nd(:)];

% Calculate Mie scattering
include_incident_field = double(true);
[efield_total, hfield_total] = mie_series(num_terms_mie, include_incident_field, ...
    refractive_index_sphere, refractive_index_background, efield_amplitude_incident, ...
    radius_sphere, wavelength, vertices_grid);
include_incident_field = double(false);
[efield_scattered, hfield_scattered] = mie_series(num_terms_mie, include_incident_field, ...
    refractive_index_sphere, refractive_index_background, efield_amplitude_incident, ...
    radius_sphere, wavelength, vertices_grid);

[efield_total_x, efield_total_y, efield_total_z] = reshape_to_grid(efield_total, x_grid_nd);
[efield_scattered_z, efield_scattered_y, efield_scattered_z] = reshape_to_grid(efield_scattered, x_grid_nd);

figure1 = figure(1);
imagesc(z_grid, x_grid, abs(efield_total_x));
xlabel('z (m)');
ylabel('x (m)');
title('|Ex| (V/m)');
colorbar;

figure2 = figure(2);
plot(x_grid, abs(efield_total_x(:, floor(num_grid_points_z/2))));
xlabel('x (m)');
ylabel('|Ex| (V/m)');
grid on;

figure3 = figure(3);
plot(z_grid, abs(efield_total_x(floor(num_grid_points_x/2), :)));
xlabel('z (m)');
ylabel('|Ex| (V/m)');
grid on;

function [output_vector_x, output_vector_y, output_vector_z] = reshape_to_grid(input_vector, x_grid_nd)
    output_vector_x = zeros(size(squeeze(x_grid_nd)));
    output_vector_y = zeros(size(squeeze(x_grid_nd)));
    output_vector_z = zeros(size(squeeze(x_grid_nd)));

    output_vector_x(:) = input_vector(:, 1);
    output_vector_y(:) = input_vector(:, 2);
    output_vector_z(:) = input_vector(:, 3);
end
