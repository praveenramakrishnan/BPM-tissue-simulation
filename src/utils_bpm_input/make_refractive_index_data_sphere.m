function [refractive_index_ndgrid] = ...
        make_refractive_index_data_sphere( ...
        filename_sphere_ensemble_params, ...
        x_grid, y_grid, z_grid, ...
        filename_sphere_data)

    input_parameters = load(filename_sphere_ensemble_params);
    refractive_index_sphere = input_parameters.refractive_index_sphere;
    refractive_index_background = input_parameters.refractive_index_background;

    sphere_data = load(filename_sphere_data);
    x_sphere_center_list = sphere_data.x_sphere_center_list;
    radius_sphere_list = sphere_data.radius_sphere_list;

    refractive_index_sphere_list = refractive_index_sphere*ones(size(radius_sphere_list));

    refractive_index_ndgrid = create_refractive_index_sphere( ...
        x_grid, y_grid, z_grid, ...
        radius_sphere_list, x_sphere_center_list, ...
        refractive_index_sphere_list, ...
        refractive_index_background);
end
