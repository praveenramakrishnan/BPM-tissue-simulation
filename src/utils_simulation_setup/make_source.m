function [efield_initial] = make_source(efield_illumination_function,...
        x_grid, y_grid, z_grid, lambda, refractive_index_background)
    z_grid_reduced = linspace(z_grid(1), 0, 2);
    [x_ndgrid, y_ndgrid, z_ndgrid] = ndgrid(x_grid, y_grid, z_grid_reduced);
    efield_illumination_function = str2func(efield_illumination_function);
    source_field = efield_illumination_function(x_ndgrid, y_ndgrid, z_ndgrid, lambda, refractive_index_background);
    efield_initial = source_field{1};
    efield_initial = efield_initial(:, :, 1);
end
