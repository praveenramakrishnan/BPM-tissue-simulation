function [efield_propagated_bpm] = solve_bpm(filename_input_parameters)
    load(filename_input_parameters);

    [refractive_index_data] = create_refractive_index_data_sphere(refractive_index_background, refractive_index_sphere, ...
        lambda, num_lambdas_lateral, length_along_z, radius_sphere, num_points_per_lambda); 
    refractive_index_average = mean(refractive_index_data(:));
    [efield_initial, x_grid, y_grid] = read_efield_input_data(filename_efield_input);

    [x_grid_nd, y_grid_nd] = ndgrid(x_grid, y_grid);

    [frequency_x_grid, frequency_y_grid] = compute_fft_frequencies(x_grid, y_grid);

    [frequency_x_grid_nd, frequency_y_grid_nd] = ndgrid(frequency_x_grid, frequency_y_grid);

    frequency_z_grid_nd = sqrt(refractive_index_average^2/lambda^2 ...
        - frequency_x_grid_nd.^2 - frequency_y_grid_nd.^2);

    num_grid_points_x = length(x_grid);
    num_grid_points_y = length(y_grid);
    num_grid_points_z_ri = size(refractive_index_data, 3);

    efield_propagated_bpm = run_beam_propagation(efield_initial, num_grid_points_x, num_grid_points_y, ...
        num_grid_points_z, step_size_z, frequency_z_grid_nd, apply_phase_correction, num_grid_points_z_ri, ...
        refractive_index_data, lambda, refractive_index_average, use_gpu_arrays);

    save(filename_variables_solve_bpm, 'efield_propagated_bpm','efield_initial', 'x_grid', 'y_grid');
end
