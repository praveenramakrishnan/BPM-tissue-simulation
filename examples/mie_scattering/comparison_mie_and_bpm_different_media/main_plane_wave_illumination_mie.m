% close all; clear all;
addpath('utils_mie_scattering');
filename_input_parameters = 'input_parameters_beam_propagation_fft_mie.mat';
load(filename_input_parameters);

if ~exist(output_directory_mie)
    mkdir(output_directory_mie);
end
  
display(append('Beginning of Mie scattering simulation'));
display(append('Output directory: ', output_directory_mie));

efield_illumination_function_name = "efield_plane_wave_mie";
efield_illumination_function = str2func(efield_illumination_function_name);

filename_output_efield = append(output_directory_mie, 'efield_illumination.mat');

if recompute_mie
    x_grid = linspace(-width_lateral/2, width_lateral/2, num_grid_points_width);
    y_grid = x_grid;
    z_grid = linspace(-length_along_z/2, length_along_z/2, 2);

     efield_illumination = efield_illumination_function(x_grid, y_grid, z_grid, lambda, ...
        refractive_index_background, refractive_index_sphere, radius_sphere, efield_amplitude_incident_field);

    save(filename_output_efield, "efield_illumination", "x_grid", "y_grid", "z_grid", "-v7.3");
else
    x_grid = load(filename_output_efield).x_grid;
    y_grid = load(filename_output_efield).y_grid;
    z_grid = load(filename_output_efield).z_grid;
    efield_illumination = load(filename_output_efield).efield_illumination;
end

if show_plots_mie
    figure1 = figure('Name', 'Mie illumination: initial field x');
    plot(x_grid, abs(squeeze(efield_illumination{1}(:, floor(num_grid_points_width/2), 1))));
    xlabel("x [m]");
    ylabel("$|E|$", "Interpreter", "latex");
    grid on;
    saveas(figure1, output_directory_mie+"figure_efield_illumination_z_0.png");

    figure2 = figure('Name', 'Mie illuminiation: final field x');
    plot(x_grid, abs(squeeze(efield_illumination{1}(:, floor(num_grid_points_width/2), 2))));
    xlabel("x [m]");
    ylabel("$|E|$", "Interpreter", "latex");
    grid on;
    saveas(figure2, append(output_directory_mie, 'figure_efield_illumination_z_minus_', ...
        num2str(length_along_z), '.png'));

    figure3 = figure('Name', 'Mie illumination: initial field 2D');
    imagesc(x_grid, y_grid, abs(squeeze(efield_illumination{1}(:, :, 1))));
    xlabel("x [m]");
    ylabel("y [m]");
    colorbar;
    saveas(figure3, append(output_directory_mie, 'figure_efield_illumination_z_0_2d.png'));

    figure4 = figure('Name', 'Mie illumination: final field 2D');
    imagesc(x_grid, y_grid, abs(squeeze(efield_illumination{1}(:, :, 2))));
    xlabel("x [m]");
    ylabel("y [m]");
    colorbar;
    saveas(figure4, append(output_directory_mie, 'figure_efield_illumination_z_minus_', ...
        num2str(length_along_z), '_2d.png'));
end
display(append('End of Mie scattering simulation'));

function [efield_illumination] = efield_plane_wave_mie(x_grid, y_grid, z_grid, lambda, ...
        refractive_index_background, refractive_index_sphere, radius_sphere, efield_amplitude_incident_field )

    size_parameter = 2*pi*radius_sphere*refractive_index_background/lambda;
    wisecombe_number = size_parameter + 4.05*size_parameter^(1/3) + 2;
    num_terms_mie = max(50, 2*wisecombe_number);

    num_grid_points_x = length(x_grid);
    num_grid_points_y = length(y_grid);

    [x_grid_nd, y_grid_nd, z_grid_nd] = ndgrid(x_grid, y_grid, z_grid);
    vertices_grid = [x_grid_nd(:), y_grid_nd(:), z_grid_nd(:)];

    [efield_total, hfield_total] = mie_series(num_terms_mie, 1, refractive_index_sphere, ...
        refractive_index_background, efield_amplitude_incident_field, radius_sphere, ... 
        lambda, vertices_grid);

    efield_illumination_x = zeros(size(squeeze(x_grid_nd)));
    efield_illumination_y = zeros(size(squeeze(x_grid_nd)));
    efield_illumination_z = zeros(size(squeeze(x_grid_nd)));

    efield_illumination_x = efield_total(:, 1);
    efield_illumination_y = efield_total(:, 2);
    efield_illumination_z = efield_total(:, 3);

    efield_illumination{1} = reshape(efield_illumination_x, num_grid_points_x, num_grid_points_y, 2);

    efield_illumination{2} = reshape(efield_illumination_y, num_grid_points_x, num_grid_points_y, 2);
end
