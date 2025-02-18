function [material_id_ndgrid] = create_material_id_spheres(x_ndgrid, y_ndgrid, z_ndgrid, radius_sphere_list, x_center_sphere_list)
    num_spheres = numel(radius_sphere_list);
    material_id_ndgrid = zeros(size(x_ndgrid));
    for iloop_sphere = 1:num_spheres
        radius_sphere = radius_sphere_list(iloop_sphere);
        x_center_sphere = x_center_sphere_list(iloop_sphere, :);

        inside_sphere = (x_ndgrid - x_center_sphere(1)).^2 ...
            + (y_ndgrid - x_center_sphere(2)).^2 ...
            + (z_ndgrid - x_center_sphere(3)).^2 ...
            <= radius_sphere^2;

        material_id_ndgrid(inside_sphere) = iloop_sphere;
    end
end


