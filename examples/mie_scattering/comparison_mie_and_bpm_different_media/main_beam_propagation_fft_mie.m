clear all; close all;
addpath('utils_mie_scattering');
addpath('./utils_bpm')

num_refractive_indices = 13;
num_radii = 4;

for iloop_ri = 1:num_refractive_indices 
    display(append('Begin ri loop number ', num2str(iloop_ri)));
    for iloop_radius=1:num_radii
        display(append('Begin radii loop number ', num2str(iloop_radius)));
        % Run the input parameters file
        filename_input_parameters = 'input_parameters_beam_propagation_fft_mie.mat';
        run(strrep(filename_input_parameters, '.mat', '.m'));

        % Run Mie scattering calculation
        run('./main_plane_wave_illumination_mie.m');

        % Create output directory if it doesn't exist
        if ~exist(output_directory_bpm)
            mkdir(output_directory_bpm);
        end

        display("Begin BPM simulation")
        display(append('Output directory: ', output_directory_bpm));

        if recompute_bpm
            tic;
            efield_propagated_bpm = solve_bpm(filename_input_parameters);
            time_toc = toc;
            display(append('Elapsed simulation time = ', num2str(time_toc), ' seconds'));
        end

        post_process_bpm_solution(filename_input_parameters, filename_variables_solve_bpm);

        display('End of BPM simulation');
    end
end
