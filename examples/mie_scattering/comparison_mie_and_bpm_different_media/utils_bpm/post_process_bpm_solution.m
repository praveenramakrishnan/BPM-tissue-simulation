function post_process_bpm_solution(filename_input_parameters, filename_variables_solve_bpm)
    load(filename_input_parameters);
    load(filename_variables_solve_bpm);

    % Load output from focusing code
    data_efield_illumination = load(filename_efield_input);
    efield_mie_input = data_efield_illumination.efield_illumination{1}(:, :, 2);

    % Zoom into a central portion of the data 
    Nx_out = num_grid_points_width;
    Ny_out = num_grid_points_width;
    efield_initial_stripped = strip_vectors(efield_initial, Nx_out, Ny_out);
    efield_mie_input_stripped = strip_vectors(efield_mie_input, Nx_out, Ny_out);
    efield_propagated_bpm_stripped = strip_vectors(efield_propagated_bpm, Nx_out, Ny_out);
    x_grid_stripped = strip_vectors(x_grid, 1, Nx_out);
    y_grid_stripped = strip_vectors(y_grid, 1, Ny_out);

    % Compute L2 error on entire grid
    relative_error_l2 = sqrt(sum(abs(efield_propagated_bpm(:) - efield_mie_input(:)).^2) ...
        /sum(abs(efield_mie_input(:)).^2));
    display(relative_error_l2);

    % Compute L2 error on reduced grid
    relative_error_l2_stripped = sqrt(sum(abs(efield_propagated_bpm_stripped(:) - ...
        efield_mie_input_stripped(:)).^2)/sum(abs(efield_mie_input(:)).^2));
    display(relative_error_l2_stripped);

    % Save data
    save(filename_efield_output, 'efield_propagated_bpm', 'efield_mie_input', 'x_grid', 'y_grid', 'relative_error_l2', 'refractive_index_sphere', 'radius_sphere');
    writematrix([length_along_z, num_lambdas_lateral, relative_error_l2], ...
        append(output_directory_bpm, 'relative_error_l2_mie_field.txt'), ...
        'Delimiter', 'tab');
    writematrix([length_along_z, num_lambdas_lateral, relative_error_l2_stripped], ...
        append(output_directory_bpm, 'relative_error_l2_mie_field_reduced_num_points_', ...
        num2str(Nx_out), 'x', num2str(Ny_out),'.txt'), 'Delimiter', 'tab');

    % Plot data on reduced grid
    if show_plots_bpm
        figure('Name', 'BPM solution: initial field x');
        plot(x_grid_stripped, abs(efield_initial_stripped(:, Ny_out/2)), 'DisplayName', 'Field at input plane along x');
        xlabel('x (m)');
        ylabel('|E| (V/m)');
        title('Magnitude of electric field at the input plane');
        legend;
        grid on;

        figure('Name', 'BPM solution: final field x');
        plot(x_grid_stripped, abs(efield_mie_input_stripped(:, Ny_out/2)), 'DisplayName', 'Mie solution');
        hold on;
        plot(x_grid_stripped, abs(efield_propagated_bpm_stripped(:, Ny_out/2)), 'DisplayName', 'Beam propagation method');
        hold off;
        xlabel('x (m)');
        ylabel('|E| (V/m)');
        title('Magnitude of electric field at the output plane');
        legend;
        grid on;

        figure('Name', 'BPM solution: initial field y');
        plot(y_grid_stripped, abs(efield_initial_stripped(Nx_out/2, :)), 'DisplayName', 'Field at input plane along y');
        xlabel('y (m)');
        ylabel('|E| (V/m)');
        title('Magnitude of electric field at the input plane');
        legend;
        grid on;

        figure('Name', 'BPM solution: final field y');
        plot(y_grid_stripped, abs(efield_mie_input_stripped(Nx_out/2, :)), 'DisplayName', 'Mie solution');
        hold on;
        plot(y_grid_stripped, abs(efield_propagated_bpm_stripped(Nx_out/2, :)), 'DisplayName', 'Beam propagation method');
        hold off;
        xlabel('y (m)');
        ylabel('|E| (V/m)');
        title('Magnitude of electric field at the output plane');
        legend;
        grid on;

        % figure(11);
        % plot(x_grid_stripped, abs(efield_propagated_bpm_stripped(:, Ny_out/2) - ...
        %     efield_mie_input_stripped(:, Ny_out/2)), 'DisplayName', 'Difference in calculated fields');
        % xlabel('x (m)');
        % ylabel('$|E_{BPM} - E_{mie}|$ (V/m)', 'Interpreter', 'Latex');
        % title('Magnitude of difference between electric fields using BPM and Debye-Wolf methods');
        % legend;
        % grid on;

        % figure(12);
        % plot(y_grid_stripped, abs(efield_propagated_bpm_stripped(Nx_out/2, :) - ...
        %     efield_mie_input_stripped(Nx_out/2, :)), 'DisplayName', 'Difference in calculated fields');
        % xlabel('y (m)');
        % ylabel('$|E_{BPM} - E_{mie}|$ (V/m)', 'Interpreter', 'Latex');
        % title('Magnitude of difference between electric fields using BPM and Debye-Wolf methods');
        % legend;
        % grid on;

        % figure(101);
        % imagesc(x_grid, y_grid, abs(efield_initial));
        % xlabel('x (m)');
        % ylabel('y (m)');
        % colorbar;
        % title("Input electric field (V/m)");

        % figure(102);
        % imagesc(x_grid_stripped, y_grid_stripped, abs(efield_mie_input_stripped));
        % xlabel('x (m)');
        % ylabel('y (m)');
        % colorbar;
        % title("Output electric field Mie solution (V/m)");

        % figure(103);
        % imagesc(x_grid_stripped, y_grid_stripped, abs(efield_propagated_bpm_stripped));
        % xlabel('x (m)');
        % ylabel('y (m)');
        % title("Output electric field BPM (V/m)");
        % colorbar;
    end
end
