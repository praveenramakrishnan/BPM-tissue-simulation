% Add path to library functions
clear all; close all;
run('../tdms_location.m');
addpath(append(tdms_root, 'tdms/tests/system/data/input_generation/matlab'));
addpath(genpath('../../../src/utils_simulation_setup/'));
addpath('../../../src/utils_pstd_solution/');

% Parameters
recompute_refractive_index = true;
save_refractive_index = true;
recompute_material_grid = true;
recompute_efield_initial_pstd = true;
recompute_pstd_setup = true;
recompute_pstd_solution = true;
recompute_post_process_pstd_solution = true;
simulate_background = false;

% Load input parameters
filename_input_parameters_pstd = 'input_file_sphere_pstd.m';
run(filename_input_parameters_pstd);

% Lengths as multiples of wavelength
num_lambda_x = num2str(length_along_x/lambda);
num_lambda_y = num2str(length_along_y/lambda);
num_lambda_z = num2str(length_along_z/lambda);

% Name of the directory to which output files are saved
output_directory = append('./simulations_pstd_', source_type, '_', ...
    num_lambda_x, 'lx', num_lambda_y, 'lx', num_lambda_z, 'l/');

if simulate_background
    output_directory = append(output_directory, 'simulation_background/');
end

% Filenames
filename_grid = append(output_directory, 'grid_data.mat');
filename_material_grid = append(output_directory, 'material_grid_data.mat');
filename_refractive_index_data = append(output_directory, 'refractive_index_data_3d.mat');
filename_efield_initial_pstd = append(output_directory, 'efield_initial_pstd.mat');
filename_pstd_setup = append(output_directory, 'setup_pstd_simulation.mat');
filename_pstd_iterations = append(output_directory, 'pstd_simulation_iterations.txt');
filename_pstd_output_data =  append(output_directory, 'pstd_output_data.mat');
filename_pstd_output_efield =  append(output_directory, 'efield_pstd.mat');

% Create output directory
if ~exist(output_directory)
    mkdir(output_directory);
end

% Move created data files to output directory
if ~isempty(dir('*.mat'))
    movefile('*.mat', output_directory);
end

% Make grid
display("Start make grid");
tic;
[x_grid, y_grid, z_grid, lambda] = make_grid_pstd(filename_input_parameters_pstd);
save(filename_grid, 'x_grid', 'y_grid', 'z_grid', 'lambda', '-v7.3');
toc;
display("End make grid");

% Parameters defining the sphere
radius_sphere_list = [radius_sphere];
x_sphere_center_list = [[0, 0, 0]];
if simulate_background
    refractive_index_sphere_list = [refractive_index_background];
else
    refractive_index_sphere_list = [refractive_index_sphere];
end

display("Start make refractive index");
tic;
if recompute_refractive_index
    refractive_index_data = create_refractive_index_sphere( ...
        x_grid, y_grid, z_grid, ...
        radius_sphere_list, x_sphere_center_list, ...
        refractive_index_sphere_list, ...
        refractive_index_background);
    if save_refractive_index
        save(filename_refractive_index_data, 'refractive_index_data', ...
            'x_grid', 'y_grid', 'z_grid', '-v7.3');
    end
elseif exist(filename_refractive_index_data)
    refractive_index_data = load(filename_refractive_index_data).refractive_index_data;
else
    error(append('File ', filename_refractive_index_data, ' not found.'));
end
toc;
display("End make refractive index");

display('Begin creating material grid file');
tic;
if recompute_material_grid
    [composition_matrix, material_matrix] = ...
        create_material_grid_sphere(filename_input_parameters_pstd, ...
        radius_sphere_list, x_sphere_center_list, ...
        refractive_index_sphere_list);
    save(filename_material_grid, 'material_matrix', 'composition_matrix', '-v7.3');
end
toc;
display('End creating material grid file');

display("Start make PSTD source");
tic;
if recompute_efield_initial_pstd
    % use iteratefdtd_matrix to set up illumination file
    iteratefdtd_matrix(filename_input_parameters_pstd, 'illsetup', ...
    filename_efield_initial_pstd, filename_material_grid, '');
end
if exist(filename_efield_initial_pstd)
    efield_initial_pstd = squeeze(load(filename_efield_initial_pstd).Ksource(1, :, :));
end
toc;
display("End make PSTD source");

display("Start make PSTD simulation setup");
tic;
if recompute_pstd_setup
    display('Begin iterate fdtd calculation for sphere file setup');
    % use iteratefdtd_matrix to set up file for tdms execution
    iteratefdtd_matrix(filename_input_parameters_pstd, 'filesetup', filename_pstd_setup,...
        filename_material_grid, filename_efield_initial_pstd);
end
toc;
display("End make PSTD simulation setup");

% Run the PSTD simulation
display("Begin PSTD solution");
if recompute_pstd_solution
    system(append('echo > ', filename_pstd_iterations));
    run_tdms_command_sphere = append('tdms ', filename_pstd_setup, ' ', ...
        filename_pstd_output_data, ' >> ', ...
        filename_pstd_iterations);
    system(run_tdms_command_sphere);
end
display("End PSTD solution");

% Clean up directory and move all outputs to right folder
if ~isempty(dir('*.mat'))
    movefile('*.mat', output_directory);
end

% Post-process the PSTD solution
display("Begin post process PSTD solution");
if recompute_post_process_pstd_solution
    data_pstd_output = load(filename_pstd_output_data);
    % Reshape the output field data
    efield_samples = data_pstd_output.campssample;
    [i, j, k] = size(ii);
    efield_propagated_pstd = zeros(i, j, k);
    efield_propagated_pstd(:) = efield_samples(:);
    % Save the solution
    save(filename_pstd_output_efield, 'efield_propagated_pstd', 'x_grid', 'y_grid');
end
display("End post process PSTD solution");

% Plot results
figure_incident = figure(1);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_initial_pstd));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
colorbar;
title("Electric field magnitude (V/m)");
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
saveas(figure_incident, append(output_directory, ...
        'figure_incident_2d.png'));

figure_pstd_2d = figure(2);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_propagated_pstd));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
title("Electric field magnitude (V/m)");
colorbar;
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
saveas(figure_pstd_2d, append(output_directory, ...
    'figure_pstd_2d.png'));

figure_pstd_1d_along_x = figure(3);
plot(1e6*x_grid, abs(efield_propagated_pstd(:, find(y_grid==0))), 'DisplayName', 'PSTD');
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*x_grid), max(1e6*x_grid)]);
legend;
grid on;
saveas(figure_pstd_1d_along_x, append(output_directory, ...
    'figure_pstd_along_x.png'));

figure_pstd_1d_along_y = figure(4);
plot(1e6*x_grid, abs(efield_propagated_pstd(find(x_grid==0), :)), 'DisplayName', 'PSTD');
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*y_grid), max(1e6*y_grid)]);
legend;
grid on;
saveas(figure_pstd_1d_along_y, append(output_directory, ...
    'figure_pstd_along_x.png'));
