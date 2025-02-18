function [composition_matrix] = create_composition_matrix_cylinders(x_ndgrid,z_ndgrid, radius_cylinder_list, x_center_cylinder_list)
   num_cylinders = numel(radius_cylinder_list);
   composition_matrix = [];
   for iloop_cylinder = 1:num_cylinders
        radius_cylinder = radius_cylinder_list(iloop_cylinder);
        x_center_cylinder = x_center_cylinder_list(iloop_cylinder, :);

        inside_cylinder = (x_ndgrid - x_center_cylinder(1)).^2 ...
            + (z_ndgrid - x_center_cylinder(3)).^2 ...
            <= radius_cylinder^2;

        [ii, kk] = ind2sub(size(inside_cylinder), find(inside_cylinder>0));
        material_id = iloop_cylinder*ones(size(ii));
        composition_matrix = [composition_matrix; ii(:), ones(size(ii(:))), kk(:), material_id];
   end
end
