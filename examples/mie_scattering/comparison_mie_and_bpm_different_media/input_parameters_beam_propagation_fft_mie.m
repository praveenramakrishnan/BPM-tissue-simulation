% Parameters
recompute_mie = true;
recompute_bpm = true;
lambda = 1300e-9;
length_along_z = 50*lambda;
num_lambdas_lateral = 150;
width_lateral = num_lambdas_lateral*lambda;
num_points_per_lambda = 10;
num_grid_points_width = num_lambdas_lateral*num_points_per_lambda;
refractive_index_background = 1.3333;
% refractive_index_sphere = 1.4;
refractive_index_sphere = (iloop_ri-1)*0.025 + 1.2;
% radius_sphere = 5*lambda;
radius_sphere = (iloop_radius*lambda)*5;
efield_amplitude_incident_field = 1; % For mie solver
apply_phase_correction = true;
show_plots_mie = false; % For mie solver
show_plots_bpm = false;
use_gpu_arrays = false;

step_size_z = lambda/num_points_per_lambda;
num_grid_points_z = length_along_z/step_size_z;

output_directory_mie = append('artefacts_plane_illumination_mie_width_', ...
    num2str(num_lambdas_lateral), 'xlambda_depth_', num2str(length_along_z/1e-6), ...
    'e-6_npoints_', num2str(num_grid_points_width), ...
    '_radius_',num2str(radius_sphere/lambda), 'xlambda_nsph_', num2str(refractive_index_sphere), ...
    '_nbg_', num2str(refractive_index_background), '/');

output_directory_bpm = append('artefacts_bpm_fft_width_', num2str(num_lambdas_lateral), ...
    'xlambda_depth_', num2str(length_along_z/1e-6), 'e-6_npoints_', ...
    num2str(num_grid_points_width), '_radius_', num2str(radius_sphere/lambda), 'xlambda_nsph_', ...
    num2str(refractive_index_sphere), '_nbg_', num2str(refractive_index_background), ...
    '_numiter_', num2str(num_grid_points_z), '/');

filename_efield_input = append(output_directory_mie, 'efield_illumination.mat');
filename_efield_output = append(output_directory_bpm, 'efield_bpm.mat');
filename_variables_solve_bpm = append(output_directory_bpm, 'variables_solve_bpm.mat');

save('input_parameters_beam_propagation_fft_mie.mat');
