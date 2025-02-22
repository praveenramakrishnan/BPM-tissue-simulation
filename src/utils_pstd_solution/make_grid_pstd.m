function [x_grid, y_grid, z_grid, lambda] = make_grid_pstd(filename_input_parameters_pstd)
    [x_grid, y_grid, z_grid, lambda] = fdtd_bounds(filename_input_parameters_pstd);

    if numel(y_grid) < 2
        y_grid = 0;
    end
end
