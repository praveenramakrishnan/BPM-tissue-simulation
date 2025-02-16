function [efield_total_x, efield_total_y, efield_total_z] = compute_mie_solution(x_grid, y_grid, z_grid, ...
        lambda, radius_sphere, efield_amplitude_incident_field, refractive_index_sphere, ...
        refractive_index_background, num_terms_mie_min) 

    % Set the number of terms in mie series based on Wiscombe number
    size_parameter = 2*pi*radius_sphere*refractive_index_background/lambda;
    wisecombe_number = size_parameter + 4.05*size_parameter^(1/3) + 2;
    num_terms_mie = max(num_terms_mie_min, 2*wisecombe_number);

    % Get grid vertices
    [x_ndgrid, y_ndgrid, z_ndgrid] = ndgrid(x_grid, y_grid, z_grid);
    vertices_grid = [x_ndgrid(:), y_ndgrid(:), z_ndgrid(:)];

    % Run Mie series solver
    include_incident_field = double(true); % true for total field
    [efield_total, hfield_total] = mie_series(num_terms_mie, include_incident_field, refractive_index_sphere, ...
        refractive_index_background, efield_amplitude_incident_field, radius_sphere, ... 
        lambda, vertices_grid);

    [efield_total_x, efield_total_y, efield_total_z] = reshape_to_grid(efield_total, x_ndgrid);
end

function [output_vector_x, output_vector_y, output_vector_z] = reshape_to_grid(input_vector, x_ndgrid)
    % The output from mie_series calculated is to be reformatted to the shape of the grid.
    output_vector_x = zeros(size(squeeze(x_ndgrid)));
    output_vector_y = zeros(size(squeeze(x_ndgrid)));
    output_vector_z = zeros(size(squeeze(x_ndgrid)));

    output_vector_x(:) = input_vector(:, 1);
    output_vector_y(:) = input_vector(:, 2);
    output_vector_z(:) = input_vector(:, 3);
end
