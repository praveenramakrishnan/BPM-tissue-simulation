function [refractive_index_data] = create_average_refractive_index_data(refractive_index_average, ...
        lambda, num_lambdas_lateral, length_along_z, num_points_per_lambda) 
    % The discretization levels matching POL refractive index data
    delta_region_x = lambda/num_points_per_lambda;
    delta_region_y = lambda/num_points_per_lambda;
    delta_region_z = lambda/num_points_per_lambda;

    width_region_x = num_lambdas_lateral*lambda;
    width_region_y = num_lambdas_lateral*lambda;

    num_points_x = floor(width_region_x/delta_region_x);
    num_points_y = floor(width_region_y/delta_region_y);
    num_points_z = floor(length_along_z/delta_region_z);

    refractive_index_data= refractive_index_average*ones(num_points_x, num_points_y, num_points_z);
end

