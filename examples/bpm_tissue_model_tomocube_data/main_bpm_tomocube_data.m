clear all; close all;
addpath(genpath('../../src/'));

% Parameters
recompute_efield_initial = true;
recompute_refractive_index = true;
save_refractive_index = true;
recompute_bpm_solution = true;

simulate_background = false;

show_plot = true;
z_grid_measurement = 0e-6; % z coordinate at which BPM output is computed

% Load BPM input parameters
filename_input_parameters = 'input_parameters_bpm';
run(filename_input_parameters);
data_input = load(filename_input_parameters);

% Lengths as multiples of wavelength
num_lambda_x = num2str(length_along_x/lambda);
num_lambda_y = num2str(length_along_y/lambda);
num_lambda_z = num2str(length_along_z/lambda);

% The portion of refractive index data used as defined in 'src/utils_simulation_setup/make_refractive_index_data_tomocube.m'
scattering_geometry = '3';

% Select source type (plane wave or focused beam)
% source_type = 'plane';
source_type = 'focused';

% Name of the directory to which output files are saved
output_directory = append('./simulations_bpm_', source_type, '_', ...
    num_lambda_x, 'lx', num_lambda_y, 'lx', num_lambda_z, 'l_', ...
    'nslices_', num2str(num_region_slices), '/');

if simulate_background
    output_directory = append(output_directory, 'simulations_background/');
end

% Filenames
filename_grid = append(output_directory, 'grid_data');
filename_refractive_index_data = append(output_directory, 'refractive_index_data_3d');
filename_efield_initial = append(output_directory, 'efield_initial');
filename_bpm_efield_output = append(output_directory, 'efield_bpm.mat');

filename_input_parameters_tomocube = 'input_parameters_tomocube_data';

% Read filename for tomocube refractive index data
run(filename_input_parameters_tomocube);
filename_tomocube_data = load(filename_input_parameters_tomocube).filename_tomocube_data;

% Select the function for source computation based of source_type
if strcmp(source_type, 'plane')
    efield_illumination_function = 'efield_plane';
else
    efield_illumination_function = 'efield_focused_gauss_rotsim_1D';
end

% BPM simulation without phase correction is same as simulating
% homogeneous region with average refractive index.
if simulate_background
    apply_phase_correction = false;
    filename_bpm_efield_output = append(output_directory, 'efield_background.mat');
else
    apply_phase_correction = true;
end

% Create output directory
if ~exist(output_directory)
    mkdir(output_directory);
end

refractive_index_background = load(filename_input_parameters_tomocube).refractive_index_background;

% Move created data files to output directory
if ~isempty(dir('*.mat'))
    movefile('*.mat', output_directory);
end

display("Begin BPM simulation");
display(append('Output directory: ', output_directory));

display("Start make grid");
[x_grid, y_grid, z_grid, lambda] = make_grid(append(output_directory, filename_input_parameters));
save(filename_grid, 'x_grid', 'y_grid', 'z_grid', 'lambda', '-v7.3');
display("End make grid");

display("Start make source");
if recompute_efield_initial
    efield_initial = make_source(efield_illumination_function, ...
        x_grid, y_grid, z_grid, lambda, refractive_index_background);
    save(filename_efield_initial, 'efield_initial', '-v7.3');
else
    efield_initial = load(filename_efield_initial).efield_initial;
end
display("End make source");

% Run BPM code for each slice
efield_initial_slice = efield_initial;
for iloop_slice = 1:num_region_slices

    display(append('Running slice loop ', num2str(iloop_slice)));

    z_grid_slice = slice_z_grid(z_grid, num_region_slices, iloop_slice);

    display("Start make refractive index");
    refractive_index_data = get_refractive_index_data( ...
        data_input, filename_refractive_index_data, ...
        scattering_geometry, refractive_index_background, ...
        x_grid, y_grid, z_grid_slice, iloop_slice, ...
        filename_tomocube_data, recompute_refractive_index, ...
        simulate_background, save_refractive_index);
    display("End make refractive index");

    length_along_z = z_grid_slice(...
        abs(z_grid_slice-z_grid_measurement) == ...
        min(abs(z_grid_slice-z_grid_measurement)) ...
        ) - z_grid_slice(1);

    display("Start compute BPM solution");
    if recompute_bpm_solution
        tic;

        efield_propagated_bpm = run_beam_propagation( ...
            efield_initial_slice, x_grid, y_grid, ...
            num_bpm_planes, length_along_z, ...
            refractive_index_data, refractive_index_background, ...
            lambda, apply_phase_correction);

        time_toc = toc;
        display(append('Elapsed simulation time = ', num2str(time_toc), ' seconds'));
        save(filename_bpm_efield_output, 'efield_propagated_bpm', 'x_grid', 'y_grid');
    elseif exist(filename_bpm_efield_output)
        efield_propagated_bpm = load(filename_bpm_efield_output).efield_propagated_bpm;
    else
        error(append('File ', filename_bpm_efield_output, ' not found.'));
    end
    display("End compute BPM solution");

    % Save electric field of the slice to file
        save(append(filename_efield_initial, '_slice_', num2str(iloop_slice)), ...
            'efield_initial_slice', '-v7.3');

    % Update electric field of the slice to file
    efield_initial_slice = efield_propagated_bpm;
end

% Plot results
figure_incident = figure(1001);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_initial));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
colorbar;
% title("Input electric field (V/m)");
title("Electric field magnitude (V/m)");
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
saveas(figure_incident, append(output_directory, 'figure_bpm_incident_2d.png'));

figure_bpm_2d_full = figure(1002);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_propagated_bpm));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
title("Electric field magnitude (V/m)");
colorbar;
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
saveas(figure_bpm_2d_full, append(output_directory, 'figure_bpm_brain_tissue_2d.png'));

figure_bpm_1d_along_x = figure(1);
plot(1e6*x_grid, abs(efield_propagated_bpm(:, find(y_grid==0))), 'DisplayName', 'BPM');
hold off;
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*x_grid), max(1e6*x_grid)]);
legend;
grid on;
saveas(figure_bpm_1d_along_x, append(output_directory, 'figure_bpm_brain_tissue_along_x.png'));

figure_bpm_1d_along_y = figure(2);
plot(1e6*y_grid, abs(efield_propagated_bpm(find(x_grid==0), :)), 'DisplayName', 'BPM');
hold off;
xlabel('y ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*y_grid), max(1e6*y_grid)]);
legend('location', 'northeast');
grid on;
saveas(figure_bpm_1d_along_y, append(output_directory, 'figure_bpm_brain_tissue_along_y.png'));

display('End of BPM simulation');

function refractive_index_data = get_refractive_index_data( ...
        data_input, filename_refractive_index_data, ...
        scattering_geometry, refractive_index_background, ...
        x_grid, y_grid, z_grid_slice, iloop_slice, ...
        filename_tomocube_data, recompute_refractive_index, ...
        simulate_background, save_refractive_index)

    filename_refractive_index_data_slice = append( ...
        filename_refractive_index_data, ...
        '_slice_', num2str(iloop_slice));
    if simulate_background
        K_slice = length(z_grid_slice);
        refractive_index_data = refractive_index_background*ones(...
            data_input.I, data_input.J, K_slice); 
    elseif recompute_refractive_index
            refractive_index_data = make_refractive_index_data_tomocube(filename_tomocube_data, ...
                x_grid, y_grid, z_grid_slice, scattering_geometry);
            if save_refractive_index
                save(filename_refractive_index_data_slice, 'refractive_index_data',...
                    'x_grid', 'y_grid', 'z_grid_slice', '-v7.3');
            end
    else
        refractive_index_data = load(filename_refractive_index_data_slice).refractive_index_data;
    end
end
