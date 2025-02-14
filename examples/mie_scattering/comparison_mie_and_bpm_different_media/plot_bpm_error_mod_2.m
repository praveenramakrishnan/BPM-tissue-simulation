clear all; close all;
addpath('./utils_bpm');
recompute = false;
num_refractive_indices = 13;
num_radii = 4;
lambda = 1300e-9;
% length_measurement_cube = 150*lambda;
if recompute
    for iloop_ri=1:num_refractive_indices
        for iloop_radius=1:num_radii
            filename_input_parameters = 'input_parameters_plot_bpm_error.mat';
            run(strrep(filename_input_parameters, '.mat', '.m'));

            load(filename_input_parameters);
            length_measurement_cube = 5*radius_sphere;
            refractive_index_list(iloop_ri) = refractive_index_sphere;
            radius_sphere_list(iloop_radius) = radius_sphere;
            size_parameter = 2*pi*radius_sphere*refractive_index_background/lambda;
            size_parameter_list(iloop_radius) = size_parameter;
            wisecomb_number_list(iloop_radius) = size_parameter + 4.05*size_parameter^(1/3) + 2;

            bpm_solver_output = load(append(output_directory_bpm, 'efield_bpm.mat'));

            [efield_bpm_stripped, efield_mie_stripped] = ...
                strip_field_vectors(bpm_solver_output, length_measurement_cube);
            intensity_bpm_stripped = abs(efield_bpm_stripped).^2;
            intensity_mie_stripped = abs(efield_mie_stripped).^2;

            normalised_rms_error_l2_list(iloop_ri, iloop_radius) = ...
                sqrt(mean(abs(efield_bpm_stripped(:) - efield_mie_stripped(:)).^2)) ...
                /max(abs(efield_mie_stripped(:)))

            normalised_rms_error_l4_list(iloop_ri, iloop_radius) = ...
                sqrt(mean(abs(intensity_bpm_stripped(:) - intensity_mie_stripped(:)).^2)) ...
                /max(abs(intensity_mie_stripped(:)))

            relative_error_l2_list(iloop_ri, iloop_radius) = ...
                sqrt(sum(abs(efield_bpm_stripped(:) - efield_mie_stripped(:)).^2)) ...
                /sqrt(sum(abs(efield_mie_stripped(:)).^2))

            % relative_error_l2_list(iloop_ri, iloop_radius) = bpm_solver_output.relative_error_l2;

            display(iloop_ri);
            display(iloop_radius);
            % display(relative_error_l2_list(iloop_ri, iloop_radius));
            save('variables_plot_bpm_error_mod.mat');
        end
    end
else
    load('variables_plot_bpm_error_mod.mat');
end

refractive_index_list = refractive_index_list(4:num_refractive_indices-4);
relative_error_l2_list = relative_error_l2_list(4:num_refractive_indices-4, :);
normalised_rms_error_l2_list = normalised_rms_error_l2_list(4:num_refractive_indices-4, :);
normalised_rms_error_l4_list = normalised_rms_error_l4_list(4:num_refractive_indices-4, :);

figure1 = figure(1);
for iloop_radius=1:num_radii
    plot(refractive_index_list, relative_error_l2_list(:, iloop_radius), '*-', ...
    'DisplayName', append('$r_s$=', num2str(radius_sphere_list(iloop_radius)/lambda),'$\lambda_0$' ));
hold on;
end

hold off;
xlabel('refractive index of sphere');
ylabel('$\eta_1$', 'interpreter', 'latex');
legend('Interpreter', 'latex', 'location', 'north');
grid on;
xticks(refractive_index_list);
xticks('manual');
xlim('tight');
saveas(figure1, 'relative_error_vs_refractive_index.png');

figure2 = figure(2);
for iloop_radius=1:num_radii
    plot(refractive_index_list, normalised_rms_error_l2_list(:, iloop_radius), '*-', ...
    'DisplayName', append('$r_s$=', num2str(radius_sphere_list(iloop_radius)/lambda),'$\lambda_0$' ));
hold on;
end

hold off;
xlabel('refractive index of sphere');
ylabel('$\eta_1$', 'interpreter', 'latex');
legend('Interpreter', 'latex', 'location', 'north');
grid on;
xticks(refractive_index_list);
xticks('manual');
xlim('tight');
saveas(figure2, 'rms_error_efield_vs_refractive_index.png');

figure3 = figure(3);
for iloop_radius=1:num_radii
    plot(refractive_index_list, normalised_rms_error_l4_list(:, iloop_radius), '*-', ...
    'DisplayName', append('$r_s$=', num2str(radius_sphere_list(iloop_radius)/lambda),'$\lambda_0$' ));
hold on;
end

hold off;
xlabel('refractive index of sphere');
ylabel('$\eta_2$', 'interpreter', 'latex');
legend('Interpreter', 'latex', 'location', 'north');
grid on;
xticks(refractive_index_list);
xticks('manual');
xlim('tight');
saveas(figure3, 'rms_error_intensity_vs_refractive_index.png');


function [efield_bpm_stripped, efield_mie_stripped] = ...
        strip_field_vectors(bpm_solver_output, length_measurement_cube)
    efield_bpm = bpm_solver_output.efield_propagated_bpm;
    efield_mie = bpm_solver_output.efield_mie_input;
    x_grid = bpm_solver_output.x_grid;
    y_grid = bpm_solver_output.y_grid;

    length_along_x = x_grid(end) - x_grid(1);
    length_along_y = y_grid(end) - y_grid(1);
    Nx = length(x_grid);
    Ny = length(y_grid);

    Nx_out = round(length_measurement_cube/length_along_x*Nx);
    Ny_out = round(length_measurement_cube/length_along_y*Ny);

    efield_bpm_stripped = strip_vectors(efield_bpm, Nx_out, Ny_out);
    efield_mie_stripped = strip_vectors(efield_mie, Nx_out, Ny_out);
end
