function [frequency_x_grid, frequency_y_grid] = compute_fft_frequencies(x_grid, y_grid)
    num_points_grid_x = length(x_grid);
    num_points_grid_y = length(y_grid);
    step_size_x = x_grid(2) - x_grid(1);
    step_size_y = y_grid(2) - y_grid(1);

    if mod(num_points_grid_x, 2) == 1
        frequency_x_grid = linspace(-0.5/step_size_x, 0.5/step_size_x, num_points_grid_x);
    else
        frequency_x_grid(1:num_points_grid_x/2+1) = linspace(-0.5/step_size_x, 0, num_points_grid_x/2+1);
        frequency_x_grid(num_points_grid_x/2+1:num_points_grid_x) = linspace(0, 0.5/step_size_x, num_points_grid_x/2);
    end

    if mod(num_points_grid_y, 2) == 1
        frequency_y_grid = linspace(-0.5/step_size_y, 0.5/step_size_y, num_points_grid_y);
    else
        frequency_y_grid(1:num_points_grid_y/2+1) = linspace(-0.5/step_size_y, 0, num_points_grid_y/2+1);
        frequency_y_grid(num_points_grid_y/2+1:num_points_grid_y) = linspace(0, 0.5/step_size_y, num_points_grid_y/2);
    end
end
