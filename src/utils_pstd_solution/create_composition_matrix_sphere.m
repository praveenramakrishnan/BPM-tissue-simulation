function [composition_matrix] = create_composition_matrix_sphere(x_ndgrid, y_ndgrid, z_ndgrid, radius_sphere_list, x_center_sphere_list)
    num_spheres = numel(radius_sphere_list);
    composition_matrix = [];
   for iloop_sphere = 1:num_spheres
        radius_sphere = radius_sphere_list(iloop_sphere);
        x_center_sphere = x_center_sphere_list(iloop_sphere, :);

        inside_sphere = (x_ndgrid - x_center_sphere(1)).^2 ...
            + (y_ndgrid - x_center_sphere(2)).^2 ...
            + (z_ndgrid - x_center_sphere(3)).^2 ...
            <= radius_sphere^2;

        [ii, jj, kk] = ind2sub(size(inside_sphere), find(inside_sphere>0));
        material_id = iloop_sphere*ones(size(ii));
        composition_matrix = [composition_matrix; ii(:), jj(:), kk(:), material_id];
   end
end
