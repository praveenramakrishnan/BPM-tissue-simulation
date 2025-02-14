function [efield_initial, x_grid, y_grid, z_grid] = read_tdms_input_data(filename_tdms_data_input)
    data_tdms_input = load(filename_tdms_data_input);

    efield_initial_raw = squeeze(interpolate_Ksource_along_x(data_tdms_input.Ksource(1,:,:)));
    % Use the factor of 0.5 below to remove the compact source adjustment
    efield_initial = 0.5*strip_pml_layers(efield_initial_raw, data_tdms_input, "xy"); 

    x_grid_raw = data_tdms_input.grid_labels.x_grid_labels; 
    y_grid_raw = data_tdms_input.grid_labels.y_grid_labels; 
    z_grid_raw = data_tdms_input.grid_labels.z_grid_labels; 

    x_grid = strip_pml_layers(x_grid_raw, data_tdms_input, "x");
    y_grid = strip_pml_layers(y_grid_raw, data_tdms_input, "y");
    z_grid = strip_pml_layers(z_grid_raw, data_tdms_input, "z");
end

function [stripped_vector] = strip_pml_layers(input_vector, data_tdms_input, strip_direction)
    num_pml_cells_x_lower = data_tdms_input.Dxl;
    num_pml_cells_x_upper = data_tdms_input.Dxu;
    num_pml_cells_y_lower = data_tdms_input.Dyl;
    num_pml_cells_y_upper = data_tdms_input.Dyu;
    num_pml_cells_z_lower = data_tdms_input.Dzl;
    num_pml_cells_z_upper = data_tdms_input.Dzu;

    interface_along_z = data_tdms_input.interface.K0(1);

    if strip_direction == "xy"
        stripped_vector = input_vector(...
            num_pml_cells_x_upper+1:end-(num_pml_cells_x_upper+1),...
            num_pml_cells_y_upper+1:end-(num_pml_cells_y_upper+1));
    elseif strip_direction =="x"
        stripped_vector = input_vector(...
            num_pml_cells_x_upper+1:end-(num_pml_cells_x_upper+1));
    elseif strip_direction =="y"
        stripped_vector = input_vector(...
            num_pml_cells_y_upper+1:end-(num_pml_cells_y_upper+1));
    elseif strip_direction =="z"
        stripped_vector = input_vector(...
            interface_along_z:end-(num_pml_cells_z_upper+1));
    else
        error(append('Unrecognized strip direction, ', strip_direction, ...
            ', allowed values are "xy", "x", "y" and "z"'));
    end
end

function [Ksource_interpolated] = interpolate_Ksource_along_x(Ksource)
    Ksource_interpolated = Ksource;
    Ksource_interpolated(1, 2:end, :) = 0.5*(Ksource(1, 1:end-1, :) + Ksource(1, 2:end, :));
end
