function [material_matrix] = create_material_matrix_sphere( ...
        radius_list, x_center_list, refractive_index_list)

    % Constant material parameters
    mu_r = 1; nu_c = 0; omega_p = 0;
    sigma_x = 0; sigma_y = 0; sigma_z = 0;
    sigma_star_x = 0; sigma_star_y = 0; sigma_star_z = 0;
    material_matrix_constants = [mu_r, nu_c, omega_p, ...
        sigma_x, sigma_y,sigma_z, ...
        sigma_star_x, sigma_star_y, sigma_star_z];

    num_materials = numel(radius_list);
    material_matrix = [];
   for iloop_material = 1:num_materials
        radius_material = radius_list(iloop_material);
        x_center_material = x_center_list(iloop_material, :);
        refractive_index = refractive_index_list(iloop_material);

        material_id = iloop_material;

        material_matrix = [material_matrix; ...
            material_id, refractive_index^2, material_matrix_constants];
   end
end
