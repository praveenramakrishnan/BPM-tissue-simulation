clear all; close all;
recompute = false;
num_refractive_indices = 13;
num_radii = 4;
if recompute
    for iloop_ri=1:num_refractive_indices
        for iloop_radius=1:num_radii
            filename_input_parameters = 'input_parameters_plot_bpm_error.mat';
            % filename_input_parameters = 'input_parameters_beam_propagation_fft_mie.mat';
            run(strrep(filename_input_parameters, '.mat', '.m'));

            load(filename_input_parameters);
            refractive_index_list(iloop_ri) = refractive_index_sphere;
            radius_sphere_list(iloop_radius) = radius_sphere;
            size_parameter = 2*pi*radius_sphere*refractive_index_background/lambda;
            size_parameter_list(iloop_radius) = size_parameter;
            wisecomb_number_list(iloop_radius) = size_parameter + 4.05*size_parameter^(1/3) + 2;

            bpm_solver_output = load(append(output_directory_bpm, 'efield_bpm.mat'));

            relative_error_l2_list(iloop_ri, iloop_radius) = bpm_solver_output.relative_error_l2;
            save('variables_plot_bpm_error.mat');
        end
    end
else
    load('variables_plot_bpm_error.mat');
end

refractive_index_list = refractive_index_list(3:num_refractive_indices-3);
relative_error_l2_list = relative_error_l2_list(3:num_refractive_indices-3, :);
figure1 = figure(1);
for iloop_radius=1:num_radii
    plot(refractive_index_list, relative_error_l2_list(:, iloop_radius), '*-', ...
    'DisplayName', append('radius=', num2str(radius_sphere_list(iloop_radius)/lambda),'$\lambda_0$' ));
hold on;
end

hold off;
xlabel('refractive index of sphere');
ylabel('relative error');
legend('Interpreter', 'latex', 'Location', 'Best');
grid on;
xticks(refractive_index_list);
xticks('manual');
xlim('tight');
saveas(figure1, 'relative_error_vs_refractive_index.png');

% figure(1)
% semilogy(radius_sphere_list, relative_error_l2_list, '-o');
% xlabel('radius of sphere (m)');
% ylabel('L2 relative error');
% grid on;

% figure(2)
% semilogy(size_parameter_list, relative_error_l2_list, '-o', 'DisplayName', ' L2 realtive error');
% hold on;
% semilogy(size_parameter_list, wisecomb_number_list, '-*', 'DisplayName', 'Wisecomb number');
% % hold on;
% % semilogy(size_parameter_list, asymptotic_number_list, '-+', 'DisplayName', 'Asymptotic number');
% hold off;
% xlabel('size parameter');
% % ylabel('l2 relative error');
% legend;
% grid on;

% figure(3)
% plot(wisecomb_number_list, relative_error_l2_list, '-o', 'DisplayName', ' L2 realtive error');
% xlabel('Wisecomb number');
% ylabel('L2 relative error');
% legend;
% grid on;
