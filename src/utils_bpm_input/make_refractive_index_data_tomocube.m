function [refractive_index_data_3d] = make_refractive_index_data_tomocube(filename_tomocube_data, ...
    x_grid, y_grid, z_grid, scattering_geometry)

    data=load(filename_tomocube_data);

    % useful range
    index_useful_start = 54;
    index_useful_end = 90;
    y_data_range = data.yvec(index_useful_end) - data.yvec(index_useful_start);
    y_data_mean = data.yvec(mean([index_useful_end, index_useful_start]));

    % simulation geometry
    % map refractive index data onto grid
    if scattering_geometry == '1'
        x_grid_offset = 300e-6;
        z_grid_offset = 0;
    elseif  scattering_geometry == '2'
        x_grid_offset = 300e-6 - 50e-6;
        z_grid_offset = 375e-6;
    elseif  scattering_geometry == '3'
        x_grid_offset = 300e-6 -62e-6;
        z_grid_offset = 325e-6;
    end

    x_grid_interpolated = x_grid - mean(x_grid) + x_grid_offset;
    y_grid_interpolated = y_grid - mean(y_grid) + y_data_mean;
    z_grid_interpolated = z_grid - z_grid(1) + z_grid_offset;

    [x_ndgrid_interpolated, y_ndgrid_interpolated, z_ndgrid_interpolated] = ndgrid(...
        x_grid_interpolated, y_grid_interpolated, z_grid_interpolated);

    refractive_index_data_3d = zeros(size(x_ndgrid_interpolated));
    [refractive_index_ndgrid_interpolated] = interp_data_nearest(data.RI, data.xvec, data.yvec, data.zvec,...
        [...
            x_ndgrid_interpolated(:),...
            sawtooth_extend_data(y_ndgrid_interpolated(:) - y_data_mean , y_data_range) + y_data_mean,...
            z_ndgrid_interpolated(:) ...
        ]);

    refractive_index_data_3d(:) = refractive_index_ndgrid_interpolated(:);
end
