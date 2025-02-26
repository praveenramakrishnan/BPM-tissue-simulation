clear all; close all;

addpath(genpath('../../../src/'));

% Parameters
recompute_mie_solution = true;
recompute_efield_initial = true;
recompute_refractive_index = true;
save_refractive_index = true;
recompute_bpm_solution = true;

% Load input parameters
filename_input_parameters = 'input_parameters';
run(filename_input_parameters);
load(filename_input_parameters);

% Lengths as multiples of wavelength
num_lambda_x = num2str(length_along_x/lambda);
num_lambda_y = num2str(length_along_y/lambda);
num_lambda_z = num2str(length_along_z/lambda);

% Refractive index and sphere
num_refractive_indices = length(refractive_index_sphere_sweep_list);
num_radius_sphere = length(radius_sphere_sweep_list);

iloop_figures = 1;
for iloop_medium = 1:num_refractive_indices
    refractive_index_sphere = refractive_index_sphere_sweep_list(iloop_medium);
    for iloop_radius = 1:num_radius_sphere
        radius_sphere = radius_sphere_sweep_list(iloop_radius);
        % Name of the directory to which output files are saved
        output_directory = append('./simulations_bpm_and_mie_comparison_', ...
            num_lambda_x, 'lx', num_lambda_y, 'lx', num_lambda_z, 'l_', ...
            'refractive_index_', num2str(refractive_index_sphere), ...
            '_radius_sphere_', num2str(radius_sphere), '/');

        % Filenames
        filename_grid = append(output_directory, 'grid_data');
        filename_refractive_index_data = append(output_directory, 'refractive_index_data_3d.mat');
        filename_efield_initial = append(output_directory, 'efield_initial.mat');
        filename_bpm_efield_output = append(output_directory, 'efield_bpm.mat');
        filename_mie_efield_output = append(output_directory, 'efield_mie.mat');

        % Create output directory
        if ~exist(output_directory)
            mkdir(output_directory);
        end

        % Copy created data files to output directory
        if ~isempty(dir('*.mat'))
            copyfile('*.mat', output_directory);
        end

        % Make grid
        display("Start make grid");
        [x_grid, y_grid, z_grid, lambda] = make_grid(append(output_directory, filename_input_parameters));
        save(filename_grid, 'x_grid', 'y_grid', 'z_grid', 'lambda', '-v7.3');
        display("End make grid");

        % Compute Mie solution
        display("Start Mie solution");
        if recompute_mie_solution
            efield_amplitude_incident_field = 1.0;
            [efield_mie_x, efield_mie_y, efield_mie_z] = compute_mie_solution( ...
                 x_grid, y_grid, z_grid, lambda, radius_sphere, ...
                 efield_amplitude_incident_field,refractive_index_sphere, ...
                 refractive_index_background, num_terms_mie);
             efield_propagated_mie = efield_mie_x;
             save(filename_mie_efield_output, 'efield_propagated_mie', 'x_grid', 'y_grid');
        elseif exist(filename_mie_efield_output)
            efield_propagated_mie = load(filename_mie_efield_output).efield_propagated_mie;
        else
            error(append('File ', filename_mie_efield_output, ' not found.'));
        end
        display("End Mie solution");

        % Compute BPM solution
        display("Start make source");
        if recompute_efield_initial
            efield_illumination_function = 'efield_plane';
            efield_initial = make_source(efield_illumination_function, x_grid, y_grid, z_grid, lambda, refractive_index_background );
            save(filename_efield_initial, 'efield_initial', '-v7.3');
        elseif exist(filename_efield_initial)
            efield_initial = load(filename_efield_initial).efield_initial;
        else
            error(append('File ', filename_efield_initial, ' not found.'));
        end
        display("End make source");

        display("Start make refractive index");
        if recompute_refractive_index
            radius_sphere_list = [radius_sphere];
            x_sphere_center_list = [[0, 0, 0]];
            refractive_index_sphere_list = [refractive_index_sphere];
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
        display("End make refractive index");

        display("Start BPM solution");
        if recompute_bpm_solution
            apply_phase_correction = true;
            efield_propagated_bpm = run_beam_propagation( ...
                efield_initial,  x_grid, y_grid, ...
                num_bpm_planes, length_along_z, ...
                refractive_index_data, refractive_index_background, ...
                lambda, apply_phase_correction);
                save(filename_bpm_efield_output, 'efield_propagated_bpm', 'x_grid', 'y_grid');
        elseif exist(filename_bpm_efield_output)
            efield_propagated_bpm = load(filename_bpm_efield_output).efield_propagated_bpm; 
        else
            error(append('File ', filename_bpm_efield_output, ' not found.'));
        end
        display("End BPM solution");

        % Compare the solutions
        figure_incident = figure(5*iloop_figures+1);
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

        figure_mie_2d = figure(5*iloop_figures+2);
        imagesc(1e6*x_grid, 1e6*y_grid, squeeze(abs(efield_propagated_mie(:, :, end))));
        xlabel('x ($\mu$m)', 'interpreter', 'latex');
        ylabel('y ($\mu$m)', 'interpreter', 'latex');
        title("Electric field magnitude (V/m)");
        colorbar;
        xticks('manual');
        x_grid_ticks = yticks;
        xticks(x_grid_ticks);
        saveas(figure_mie_2d, append(output_directory, ...
            'figure_mie_2d.png'));

        figure_bpm_2d = figure(5*iloop_figures+3);
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

        figure_comparison_1d_along_x = figure(5*iloop_figures+4);
        plot(1e6*x_grid, abs(efield_propagated_bpm(:, find(y_grid==0))), 'DisplayName', 'BPM');
        hold on;
        plot(1e6*x_grid, abs(efield_propagated_mie(:, find(y_grid==0), end)), 'DisplayName', 'Mie');
        hold off;
        xlabel('x ($\mu$m)', 'interpreter', 'latex');
        ylabel('|E| (V/m)');
        xlim([min(1e6*x_grid), max(1e6*x_grid)]);
        legend;
        grid on;
        saveas(figure_comparison_1d_along_x, append(output_directory, ...
            'figure_comparison_along_x.png'));

        figure_bpm_1d_along_y = figure(5*iloop_figures+5);
        plot(1e6*y_grid, abs(efield_propagated_bpm(find(x_grid==0), :)), 'DisplayName', 'BPM');
        hold on;
        plot(1e6*y_grid, abs(efield_propagated_mie(find(x_grid==0), :, end)), 'DisplayName', 'Mie');
        xlabel('y ($\mu$m)', 'interpreter', 'latex');
        ylabel('|E| (V/m)');
        xlim([min(1e6*y_grid), max(1e6*y_grid)]);
        legend('location', 'northeast');
        grid on;
        saveas(figure_bpm_1d_along_y, append(output_directory, ...
            'figure_comparison_along_y.png'));
        iloop_figures = iloop_figures+1;
    end
end

% Clean the directory
if ~isempty(dir('*.mat'))
    delete('*.mat');
end
