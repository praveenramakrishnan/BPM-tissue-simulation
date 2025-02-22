clear all; close all;
addpath(genpath('../../../src/'));
run('../tdms_location.m');
addpath('../pstd_source_functions');
addpath(append(tdms_root, 'tdms/tests/system/data/input_generation/matlab'));

% Parameters
recompute_refractive_index = true;
save_refractive_index = true;
recompute_material_grid = true;
recompute_efield_initial_pstd = true;
recompute_pstd_setup = true;
recompute_pstd_solution = true;

recompute_efield_initial_bpm = true;
recompute_bpm_solution = true;

simulate_background = false;

% Load input parameters
filename_input_parameters_pstd = 'input_file_sphere_pstd.m';
run(filename_input_parameters_pstd);

% Lengths as multiples of wavelength
num_lambda_x = num2str(length_along_x/lambda);
num_lambda_y = num2str(length_along_y/lambda);
num_lambda_z = num2str(length_along_z/lambda);

% Name of the directory to which output files are saved
output_directory = append('./simulations_comparison_', source_type, '_', ...
    num_lambda_x, 'lx', num_lambda_y, 'lx', num_lambda_z, 'l/');

if simulate_background
    output_directory = append(output_directory, 'simulations_background/');
end

% Filenames
filename_grid = append(output_directory, 'grid_data.mat');
filename_pstd_material_grid = append(output_directory, 'material_grid_data.mat');
filename_refractive_index_data = append(output_directory, 'refractive_index_data.mat');
filename_pstd_efield_initial = append(output_directory, 'efield_initial_pstd.mat');
filename_pstd_setup = append(output_directory, 'setup_pstd_simulation.mat');
filename_pstd_iterations = append(output_directory, 'pstd_simulation_iterations.txt');
filename_pstd_output_data =  append(output_directory, 'pstd_output_data.mat');
filename_pstd_output_efield =  append(output_directory, 'efield_pstd.mat');

filename_bpm_efield_initial = append(output_directory, 'efield_initial.mat');
filename_bpm_efield_output = append(output_directory, 'efield_propagated_bpm.mat');

if source_type == "plane"
    efield_illumination_function = 'efield_plane';
elseif source_type == "focused"
    efield_illumination_function = 'efield_focused_gauss_rotsim_1D';
else
    error(append('Unidentified source_type ', source_type));
end
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
x_sphere_center_list = [[0, 0, z_grid(ceil(end/2))]];
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
    save(filename_pstd_material_grid, 'material_matrix', 'composition_matrix', '-v7.3');
end
toc;
display('End creating material grid file');

display("Start make PSTD source");
tic;
if recompute_efield_initial_pstd
    % use iteratefdtd_matrix to set up illumination file
    iteratefdtd_matrix(filename_input_parameters_pstd, 'illsetup', ...
    filename_pstd_efield_initial, filename_pstd_material_grid, '');
end
toc;
display("End make PSTD source");

% Save initial electric field on the input plane
% if exist(filename_pstd_efield_initial)
%     efield_initial = squeeze(load(filename_pstd_efield_initial).Ksource(1, :, :));
%     save(filename_efield_initial, 'efield_initial', '-v7.3');
% end

display("Start make PSTD simulation setup");
tic;
if recompute_pstd_setup
    display('Begin iterate fdtd calculation for sphere file setup');
    % use iteratefdtd_matrix to set up file for tdms execution
    iteratefdtd_matrix(filename_input_parameters_pstd, 'filesetup', filename_pstd_setup,...
        filename_pstd_material_grid, filename_pstd_efield_initial);
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
toc;
display("End PSTD solution");

% Post-process the PSTD solution
display("Begin post process PSTD solution");
data_pstd_output = load(filename_pstd_output_data);
% Reshape the output field data
efield_samples = data_pstd_output.campssample;
efield_propagated_pstd = zeros(num_samples_x, num_samples_y, num_samples_z);
efield_propagated_pstd(:) = efield_samples(:);
% Save the solution
save(filename_pstd_output_efield, 'efield_propagated_pstd', 'x_grid', 'y_grid');
display("End post process PSTD solution");

% BPM solution
% Start from the location where source is introduced
z_grid_bpm = z_grid(source_interface_location_z:end);
refractive_index_data_bpm = refractive_index_data(:, :, source_interface_location_z:end);
display("Start make BPM source");
if recompute_efield_initial_bpm
    efield_initial = make_source(efield_illumination_function, x_grid, y_grid, z_grid_bpm, lambda, refractive_index_background);
    save(filename_bpm_efield_initial, 'efield_initial', '-v7.3');
elseif exist(filename_bpm_efield_initial)
    efield_initial = load(filename_bpm_efield_initial).efield_initial;
else
    error(append('File ', filename_bpm_efield_initial, ' not found.'));
end
display("End make BPM source");

display("Start BPM solution");
tic;
if recompute_bpm_solution
    apply_phase_correction = true;
    % A shift of delta.z/2 to compensate for compact source condition in PSTD
    length_along_z_bpm = z_grid_bpm(end) - z_grid_bpm(1) - delta.z/2;
    efield_propagated_bpm = run_beam_propagation( ...
        efield_initial,  x_grid, y_grid, ...
        num_bpm_planes, length_along_z_bpm, ...
        refractive_index_data_bpm, refractive_index_background, ...
        lambda, apply_phase_correction);
        save(filename_bpm_efield_output, 'efield_propagated_bpm', 'x_grid', 'y_grid');
elseif exist(filename_bpm_efield_output)
    efield_propagated_bpm = load(filename_bpm_efield_output).efield_propagated_bpm; 
else
    error(append('File ', filename_bpm_efield_output, ' not found.'));
end
toc;
display("End BPM solution");

% Clean up directory and move all outputs to right folder
if ~isempty(dir('*.mat'))
    movefile('*.mat', output_directory);
end

% Compare the solutions
% Error measure
error_nrms_complex_amplitude = compute_error_normalised_rms( ...
    efield_propagated_bpm, efield_propagated_pstd);
display(error_nrms_complex_amplitude);

error_nrms_intensity = compute_error_normalised_rms( ...
    abs(efield_propagated_bpm).^2, abs(efield_propagated_pstd).^2);
display(error_nrms_intensity);

% Plot results
figure_incident = figure(1);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_initial));
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

figure_bpm_2d = figure(3);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_propagated_bpm));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
title("Electric field magnitude (V/m)");
colorbar;
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
saveas(figure_bpm_2d, append(output_directory, ...
    'figure_bpm_2d.png'));

figure_comparison_1d_along_x = figure(4);
plot(1e6*x_grid, abs(efield_propagated_bpm(:, find(y_grid==0))), 'DisplayName', 'BPM');
hold on;
plot(1e6*x_grid, abs(efield_propagated_pstd(:, find(y_grid==0))), 'DisplayName', 'PSTD');
hold off;
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*x_grid), max(1e6*x_grid)]);
legend;
grid on;
saveas(figure_comparison_1d_along_x, append(output_directory, ...
    'figure_comparison_along_x.png'));

figure_comparison_1d_along_y = figure(5);
plot(1e6*x_grid, abs(efield_propagated_bpm(find(x_grid==0), :)), 'DisplayName', 'BPM');
hold on;
plot(1e6*x_grid, abs(efield_propagated_pstd(find(x_grid==0), :)), 'DisplayName', 'PSTD');
hold off;
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*y_grid), max(1e6*y_grid)]);
legend;
grid on;
saveas(figure_comparison_1d_along_y, append(output_directory, ...
    'figure_comparison_along_x.png'));

function error_normalised_rms = compute_error_normalised_rms(vector1, vector2)
    num_size_vector = length(vector2(:));
    error_normalised_rms = sqrt(sum(abs(vector1(:) - vector2(:)).^2)/num_size_vector)/max(abs(vector2(:)));
end
