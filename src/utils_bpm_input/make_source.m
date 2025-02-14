function [efield_initial] = make_source(efield_illumination_function,...
        x_grid, y_grid, z_grid, lambda)
    z_grid_reduced = linspace(z_grid(1), 0, 2);
    [x_grid_nd, y_grid_nd, z_grid_nd] = ndgrid(x_grid, y_grid, z_grid_reduced);
    efield_illumination_function = str2func(efield_illumination_function);
    source_field = efield_illumination_function(x_grid_nd, y_grid_nd, z_grid_nd, lambda);
    efield_initial = source_field{1};
    efield_initial = efield_initial(:, :, 1);
end
