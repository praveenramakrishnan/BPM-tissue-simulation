function run_sphere_code_ensemble(x_grid, y_grid, z_grid,...
        filename_sphere_ensemble_params, ...
        filename_sphere_data)

    if ~exist(append(filename_sphere_data, '.mat')) 
        input_parameters = load(filename_sphere_ensemble_params);
        refractive_index_sphere = input_parameters.refractive_index_sphere;
        refractive_index_background = input_parameters.refractive_index_background;
        sphere_concentration = input_parameters.sphere_concentration;
        radius_sphere = input_parameters.radius_sphere;

        volume_of_domain = ...
            (x_grid(end)-x_grid(1))*(y_grid(end)-y_grid(1))*(z_grid(end)-z_grid(1));
        num_spheres = round(sphere_concentration*volume_of_domain);

        [x_sphere_center, radius_sphere] = generate_ensemble(...
            x_grid(1), x_grid(end), y_grid(1), y_grid(end), z_grid(1), z_grid(end), ...
            radius_sphere, radius_sphere, num_spheres);
        x_sphere_center_list = x_sphere_center';
        radius_sphere_list = radius_sphere';
        save(filename_sphere_data, 'x_sphere_center_list', 'radius_sphere_list');
    end
end
