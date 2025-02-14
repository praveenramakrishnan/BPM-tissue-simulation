function [refractive_index_ndgrid] = ...
        create_refractive_index_sphere(...
        x_grid, y_grid, z_grid, ...
        radius_list, x_center_list, ...
        refractive_index_sphere_list, refractive_index_background)
    
    if numel(y_grid) < 2 % For 2D grid
        [x_ndgrid, z_ndgrid] = ndgrid(x_grid, z_grid);

        refractive_index_ndgrid = create_refractive_index_sphere_2d(...
            x_ndgrid, z_ndgrid, ...
            radius_list, x_center_list, ...
            refractive_index_sphere_list, refractive_index_background);
    else
        [x_ndgrid, y_ndgrid, z_ndgrid] = ndgrid(x_grid, y_grid, z_grid);

        refractive_index_ndgrid = create_refractive_index_sphere_3d(...
            x_ndgrid, y_ndgrid, z_ndgrid, ...
            radius_list, x_center_list, ...
            refractive_index_sphere_list, refractive_index_background);
    end
end

function [refractive_index_ndgrid] = create_refractive_index_sphere_2d( ...
        x_ndgrid, z_ndgrid, ...
        radius_sphere_list, x_center_sphere_list, ...
        refractive_index_sphere_list, refractive_index_background)

    num_spheres = numel(radius_sphere_list);
    refractive_index_ndgrid = refractive_index_background*ones(size(x_ndgrid));
    for iloop_sphere = 1:num_spheres
        radius_sphere = radius_sphere_list(iloop_sphere);
        x_center_sphere = x_center_sphere_list(iloop_sphere, :);
        refractive_index = refractive_index_sphere_list(iloop_sphere);

        inside_sphere = (x_ndgrid - x_center_sphere(1)).^2 ...
            + (z_ndgrid - x_center_sphere(3)).^2 ...
            <= radius_sphere^2;

        refractive_index_ndgrid(inside_sphere) = refractive_index_sphere_list(iloop_sphere);
    end
end

function [refractive_index_ndgrid] = create_refractive_index_sphere_3d( ...
        x_ndgrid, y_ndgrid, z_ndgrid, ...
        radius_sphere_list, x_center_sphere_list, ...
        refractive_index_sphere_list, refractive_index_background)

    num_spheres = numel(radius_sphere_list);
    refractive_index_ndgrid = refractive_index_background*ones(size(x_ndgrid));

    if all(radius_sphere_list(:) == radius_sphere_list(1)) & ...
            all(refractive_index_sphere_list(:) == refractive_index_sphere_list(1))
        refractive_index_ndgrid = create_refractive_index_sphere_3d_fast(...
            x_ndgrid, y_ndgrid, z_ndgrid, ...
            radius_sphere_list, x_center_sphere_list, ...
            refractive_index_sphere_list, refractive_index_background);
    else
        refractive_index_ndgrid = create_refractive_index_sphere_3d_slow(...
            x_ndgrid, y_ndgrid, z_ndgrid, ...
            radius_list, x_center_list, ...
            refractive_index_sphere_list, refractive_index_background);
    end
end


function [refractive_index_ndgrid] = create_refractive_index_sphere_3d_slow( ...
        x_ndgrid, y_ndgrid, z_ndgrid, ...
        radius_sphere_list, x_center_sphere_list, ...
        refractive_index_sphere_list, refractive_index_background)

    num_spheres = numel(radius_sphere_list);
    refractive_index_ndgrid = refractive_index_background*ones(size(x_ndgrid));
    for iloop_sphere = 1:num_spheres
        radius_sphere = radius_sphere_list(iloop_sphere);
        x_center_sphere = x_center_sphere_list(iloop_sphere, :);
        refractive_index = refractive_index_sphere_list(iloop_sphere);

        inside_sphere = (x_ndgrid - x_center_sphere(1)).^2 ...
            + (y_ndgrid - x_center_sphere(2)).^2 ...
            + (z_ndgrid - x_center_sphere(3)).^2 ...
            <= radius_sphere^2;

        refractive_index_ndgrid(inside_sphere) = refractive_index_sphere_list(iloop_sphere);
    end
end

function [refractive_index_ndgrid] = create_refractive_index_sphere_3d_fast( ...
        x_ndgrid, y_ndgrid, z_ndgrid, ...
        radius_sphere_list, x_center_sphere_list, ...
        refractive_index_sphere_list, refractive_index_background)

    radius_sphere = radius_sphere_list(1);

    delta_x = x_ndgrid(2, 1, 1) - x_ndgrid(1, 1, 1);
    delta_y = y_ndgrid(1, 2, 1) - y_ndgrid(1, 1, 1);
    delta_z = z_ndgrid(1, 1, 2) - z_ndgrid(1, 1, 1);

    inside_sphere_local = define_inside_sphere_local_fast( ...
        delta_x, delta_y, delta_z, radius_sphere);

    inside_sphere = define_inside_sphere_fast( ...
            x_ndgrid, y_ndgrid, z_ndgrid, ...
            x_center_sphere_list, radius_sphere, ...
            inside_sphere_local);

    refractive_index_ndgrid = refractive_index_background*ones(size(x_ndgrid));
    refractive_index_ndgrid(inside_sphere>0) = refractive_index_sphere_list(1);
end

function [inside_sphere_local] = define_inside_sphere_local_fast( ...
        delta_x, delta_y, delta_z, radius_sphere)

    num_points_local_x = 2*round(radius_sphere/delta_x)+1;
    num_points_local_y = 2*round(radius_sphere/delta_y)+1;
    num_points_local_z = 2*round(radius_sphere/delta_z)+1;

    length_local_x = (num_points_local_x-1)*delta_x;
    length_local_y = (num_points_local_y-1)*delta_y;
    length_local_z = (num_points_local_z-1)*delta_z;

    x_grid_local = linspace(-0.5*length_local_x, 0.5*length_local_x, num_points_local_x);
    y_grid_local = linspace(-0.5*length_local_y, 0.5*length_local_y, num_points_local_y);
    z_grid_local = linspace(-0.5*length_local_z, 0.5*length_local_z, num_points_local_z);

    [x_ndgrid_local, y_ndgrid_local, z_ndgrid_local] = ...
        ndgrid(x_grid_local, y_grid_local, z_grid_local);

    inside_sphere_local = zeros(size(x_ndgrid_local));

    inside_sphere_local(find(...
        x_ndgrid_local.^2 + y_ndgrid_local.^2 + z_ndgrid_local.^2 ...
        < radius_sphere.^2)) = 1;
end

function [inside_sphere] = define_inside_sphere_fast( ...
        x_ndgrid, y_ndgrid, z_ndgrid, ...
        x_center_sphere_list, radius_sphere, ...
        inside_sphere_local)

    num_spheres = size(x_center_sphere_list, 1);
    x_grid = x_ndgrid(:, 1, 1);
    y_grid = y_ndgrid(1, :, 1);
    z_grid = z_ndgrid(1, 1, :);
    delta_x = x_grid(2) - x_grid(1);
    delta_y = y_grid(2) - y_grid(1);
    delta_z = z_grid(2) - z_grid(1);

    index_x_center_sphere_list = zeros(size(x_center_sphere_list));

    index_x_center_sphere_list(:, 1) = round(( ...
        x_center_sphere_list(:, 1) - min(x_grid))/delta_x)+1;
    index_x_center_sphere_list(:, 2) = round(( ...
        x_center_sphere_list(:, 2) - min(y_grid))/delta_y)+1;
    index_x_center_sphere_list(:, 3) = round(( ...
        x_center_sphere_list(:, 3) - min(z_grid))/delta_z)+1;

    index_delta_radius_x = round(radius_sphere/delta_x);
    index_delta_radius_y = round(radius_sphere/delta_y);
    index_delta_radius_z = round(radius_sphere/delta_z);

    index_x_start = index_x_center_sphere_list(:, 1) - index_delta_radius_x;
    index_x_end = index_x_center_sphere_list(:, 1) + index_delta_radius_x;
    index_y_start = index_x_center_sphere_list(:, 2) - index_delta_radius_y;
    index_y_end = index_x_center_sphere_list(:, 2) + index_delta_radius_y;
    index_z_start = index_x_center_sphere_list(:, 3) - index_delta_radius_z;
    index_z_end = index_x_center_sphere_list(:, 3) + index_delta_radius_z;

    inside_sphere = zeros(size(x_ndgrid));
    for iloop_sphere=1:num_spheres
        z_sphere_center = x_center_sphere_list(iloop_sphere, 3);
        z_sphere_min = z_sphere_center -  radius_sphere;
        z_sphere_max = z_sphere_center +  radius_sphere;
        if z_sphere_min > z_grid(1) & z_sphere_max < z_grid(end)
            inside_sphere( ...
                    (index_x_start(iloop_sphere):index_x_end(iloop_sphere)), ...
                    (index_y_start(iloop_sphere):index_y_end(iloop_sphere)), ...
                    (index_z_start(iloop_sphere):index_z_end(iloop_sphere))) ...
                = inside_sphere_local;
        elseif (z_sphere_min <= z_grid(end)) & (z_sphere_max >= z_grid(1)) 
            % Special case where the sphere is partially inside the region
            num_cells_sphere = size(inside_sphere_local, 3); 
            num_cells_left = ceil((z_grid(1)-z_sphere_min)/delta_z);
            num_cells_right = ceil((z_sphere_max-z_grid(end))/delta_z);
            index_z1 = max(0, num_cells_left)+1;
            index_z2 = num_cells_sphere - max(0, num_cells_right);
            inside_sphere( ...
                    (index_x_start(iloop_sphere):index_x_end(iloop_sphere)), ...
                    (index_y_start(iloop_sphere):index_y_end(iloop_sphere)), ...
                    (index_z_start(iloop_sphere)+index_z1-1: ...
                    index_z_end(iloop_sphere)+index_z2-num_cells_sphere)) ...
                = inside_sphere_local(:, :, index_z1:index_z2);
        end
    end
end
