function [refractive_index_data] = ...
        create_refractive_index_data_sphere(refractive_index_background, refractive_index_sphere, lambda, num_lambdas_lateral, ...
        length_along_z, radius_sphere, num_points_per_lambda) 
    % The discretization levels matching POL refractive index data
    delta_region_x = lambda/num_points_per_lambda;
    delta_region_y = lambda/num_points_per_lambda;
    delta_region_z = lambda/num_points_per_lambda;

    width_region_x = num_lambdas_lateral*lambda;
    width_region_y = num_lambdas_lateral*lambda;

    num_points_x = floor(width_region_x/delta_region_x);
    num_points_y = floor(width_region_y/delta_region_y);
    num_points_z = floor(length_along_z/delta_region_z);

    x_grid_refractive_index = linspace(-width_region_x/2, width_region_x/2, num_points_x);
    y_grid_refractive_index = linspace(-width_region_y/2, width_region_y/2, num_points_y);
    z_grid_refractive_index = linspace(-length_along_z/2, length_along_z/2, num_points_z);

    refractive_index_data= refractive_index_background*ones(num_points_x, num_points_y, num_points_z);

    for iloop=1:num_points_x
        for jloop=1:num_points_y
            for kloop=1:num_points_z
                if x_grid_refractive_index(iloop)^2 + y_grid_refractive_index(jloop)^2 ...
                        + z_grid_refractive_index(kloop)^2 < radius_sphere^2
                    refractive_index_data(iloop, jloop, kloop) = refractive_index_sphere;
                end
            end
        end
    end

end
