clear all; close all;
addpath(genpath('../../src/'));

% Parameters
recompute_efield_initial = false;
recompute_refractive_index = true;
save_refractive_index = true;
simulate_background = false;
use_refractive_index_default = false;
recompute_bpm_solution = true;

% apply_phase_correction = true;
show_plot = true;
z_grid_measurement = 0e-6;
num_grid_points_z = 150;
if use_refractive_index_default
    refractive_index_default = 1.3333;
end

filename_input_parameters = 'input_parameters_bpm';
run(filename_input_parameters);
data_input = load(filename_input_parameters);
num_lambda_x = num2str(data_input.length_along_x/data_input.lambda);
num_lambda_y = num2str(data_input.length_along_y/data_input.lambda);
num_lambda_z = num2str(data_input.length_along_z/data_input.lambda);

% source_type = 'plane';
source_type = 'focused';

output_directory = append('./simulations_bpm_', source_type, '_', ...
    num_lambda_x, 'lx', num_lambda_y, 'lx', num_lambda_z, 'l_', ...
    'nslices_', num2str(num_region_slices), '/');

if strcmp(source_type, 'plane')
    efield_illumination_function = 'efield_plane';
else
    efield_illumination_function = 'efield_focused_gauss_rotsim_1D';
end
filename_grid = append(output_directory, 'grid_data');
filename_refractive_index_data = append(output_directory, 'refractive_index_data_3d');
filename_efield_initial = append(output_directory, 'efield_initial');
filename_bpm_efield_output = append(output_directory, 'efield_bpm.mat');

filename_sphere_ensemble_params = 'input_parameters_sphere_refractive_index';

output_sphere_data_directory = 'output_data_sphere_code/';
filename_sphere_data = append(output_sphere_data_directory, ...
    'sphere_data_', num_lambda_x, 'xl_', num_lambda_y, 'xl_', num_lambda_z, 'l');

if simulate_background
    apply_phase_correction = false;
    filename_bpm_efield_output = append(output_directory, 'efield_background.mat');
else
    apply_phase_correction = true;
end

if ~exist(output_directory)
    mkdir(output_directory);
end

if ~exist(output_sphere_data_directory)
    mkdir(output_sphere_data_directory);
end

run(filename_sphere_ensemble_params);
filename_sphere_ensemble_params = append(output_directory, ...
    filename_sphere_ensemble_params);

if ~isempty(dir('*.mat'))
    system(append('mv *.mat ', output_directory));
end

display("Begin BPM simulation");
display(append('Output directory: ', output_directory));


display("Start make grid");
[x_grid, y_grid, z_grid, lambda] = make_grid(append(output_directory, filename_input_parameters));
save(filename_grid, 'x_grid', 'y_grid', 'z_grid', 'lambda', '-v7.3');
display("End make grid");

display("Start make source");
if recompute_efield_initial
    efield_initial = make_source(efield_illumination_function, x_grid, y_grid, z_grid, lambda);
    save(filename_efield_initial, 'efield_initial', '-v7.3');
else
    efield_initial = load(filename_efield_initial).efield_initial;
end
display("End make source");

if recompute_refractive_index & ~exist(append(filename_sphere_data, '.mat'))
    run_sphere_code_ensemble(x_grid, y_grid, z_grid, ...
        filename_sphere_ensemble_params, ...
        filename_sphere_data);
end

efield_initial_slice = efield_initial;
for iloop_slice = 1:num_region_slices

    display(append('Running slice loop ', num2str(iloop_slice)));

    z_grid_slice = slice_z_grid(z_grid, num_region_slices, iloop_slice);

    display("Start make refractive index");
    filename_refractive_index_data_slice = append(filename_refractive_index_data, ...
        '_slice_', num2str(iloop_slice));
    if use_refractive_index_default & simulate_background
        K_slice = length(z_grid_slice);
        refractive_index_data = refractive_index_default*ones(...
            data_input.I, data_input.J, K_slice); 
    elseif recompute_refractive_index
            refractive_index_data = make_refractive_index_data_sphere( ...
                filename_sphere_ensemble_params, ...
                x_grid, y_grid, z_grid_slice, ...
                filename_sphere_data);
            if save_refractive_index
                save(filename_refractive_index_data_slice, 'refractive_index_data',...
                    'x_grid', 'y_grid', 'z_grid_slice', '-v7.3');
            end
    else
        refractive_index_data = load(filename_refractive_index_data_slice).refractive_index_data;
    end
    display("End make refractive index");

    refractive_index_average = mean(refractive_index_data(:));

    % num_grid_points_z = find( ...
    %     abs(z_grid_slice-z_grid_measurement) == min(abs(z_grid_slice-z_grid_measurement)));

    deltaz = z_grid_slice(2) - z_grid_slice(1);
    length_along_z = z_grid_slice(...
        abs(z_grid_slice-z_grid_measurement) == ...
        min(abs(z_grid_slice-z_grid_measurement)) ...
        ) - z_grid_slice(1);
    step_size_z = length_along_z/num_grid_points_z;

    display("Start compute BPM solution");
    if recompute_bpm_solution
        tic;

        [x_grid_nd, y_grid_nd] = ndgrid(x_grid, y_grid);

        [frequency_x_grid, frequency_y_grid] = compute_fft_frequencies(x_grid, y_grid);

        [frequency_x_grid_nd, frequency_y_grid_nd] = ndgrid(frequency_x_grid, frequency_y_grid);

        frequency_z_grid_slice_nd = sqrt(refractive_index_average^2/lambda^2 ...
            - frequency_x_grid_nd.^2 - frequency_y_grid_nd.^2);

        num_grid_points_x = length(x_grid);
        num_grid_points_y = length(y_grid);
        num_grid_points_z_ri = size(refractive_index_data, 3);

        efield_propagated_bpm = run_beam_propagation(efield_initial_slice, num_grid_points_x, num_grid_points_y, ...
            num_grid_points_z, step_size_z, frequency_z_grid_slice_nd, apply_phase_correction, num_grid_points_z_ri, ...
            refractive_index_data, lambda, refractive_index_average);
        time_toc = toc;
        display(append('Elapsed simulation time = ', num2str(time_toc), ' seconds'));
        save(filename_bpm_efield_output, 'efield_propagated_bpm', 'x_grid', 'y_grid');
    else
        data_efield_bpm = load(filename_bpm_efield_output); 
        efield_propagated_bpm = data_efield_bpm.efield_propagated_bpm; 
    end
    display("End compute BPM solution");
    if simulate_background
      save(append(filename_efield_initial, ...
          '_background_slice_', num2str(iloop_slice)), ...
            'efield_initial_slice', '-v7.3');
    else
        save(append(filename_efield_initial, ...
            '_slice_', num2str(iloop_slice)), ...
            'efield_initial_slice', '-v7.3');
    end
    efield_initial_slice = efield_propagated_bpm;
end

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
if simulate_background
    saveas(figure_incident, append(output_directory, ...
        'figure_bpm_incident_2d_background.png'));
else
    saveas(figure_incident, append(output_directory, ...
        'figure_bpm_incident_2d.png'));
end

figure_bpm_2d_full = figure(1002);
imagesc(1e6*x_grid, 1e6*y_grid, abs(efield_propagated_bpm));
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('y ($\mu$m)', 'interpreter', 'latex');
title("Electric field magnitude (V/m)");
colorbar;
xticks('manual');
x_grid_ticks = yticks;
xticks(x_grid_ticks);
if simulate_background
    saveas(figure_bpm_2d_full, append(output_directory, ...
        'figure_bpm_sphere_ensemble_2d_backgorund.png'));
else 
    saveas(figure_bpm_2d_full, append(output_directory, ...
        'figure_bpm_sphere_ensemble_2d.png'));
end

figure_bpm_1d_along_x = figure(1);
plot(1e6*x_grid, abs(efield_propagated_bpm(:, find(y_grid==0))), 'DisplayName', 'BPM');
hold off;
xlabel('x ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*x_grid), max(1e6*x_grid)]);
% title(append('Magnitude of electric field at the plane z=', ...
%     num2str(z_grid_measurement), 'm'));
% xticks('manual');
% xticks(x_grid_ticks);
legend;
grid on;
if simulate_background
    saveas(figure_bpm_1d_along_x, append(output_directory, ...
        'figure_bpm_sphere_ensemble_along_x_background.png'));
    else
    saveas(figure_bpm_1d_along_x, append(output_directory, ...
        'figure_bpm_sphere_ensemble_along_x.png'));
end

figure_bpm_1d_along_y = figure(2);
plot(1e6*y_grid, abs(efield_propagated_bpm(find(x_grid==0), :)), 'DisplayName', 'BPM');
hold off;
xlabel('y ($\mu$m)', 'interpreter', 'latex');
ylabel('|E| (V/m)');
xlim([min(1e6*x_grid), max(1e6*x_grid)]);
% title(append('Magnitude of electric field at the plane z=', ...
%     num2str(z_grid_measurement), 'm'));
% xticks('manual');
% xticks(x_grid_ticks);
legend('location', 'northeast');
grid on;
if simulate_background
    saveas(figure_bpm_1d_along_y, append(output_directory, ...
        'figure_bpm_sphere_ensemble_along_y_background.png'));
else
    saveas(figure_bpm_1d_along_y, append(output_directory, ...
        'figure_bpm_sphere_ensemble_along_y.png'));
end
% Clean up directory and move all outputs to right folder
if ~isempty(dir('*.mat'))
    system(append('mv *.mat ', output_directory));
end
% if ~isempty(dir('nohup_*'))
%     system(append('mv nohup_* ', output_directory));
% end

display('End of BPM simulation');
