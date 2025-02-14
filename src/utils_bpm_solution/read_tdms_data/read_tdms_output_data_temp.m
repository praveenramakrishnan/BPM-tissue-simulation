function [efield_tdms_output] = read_tdms_output_data(filename_tdms_data_output, ...
        filename_tdms_data_input, length_measurement_cube)

    data_tdms_output = load(filename_tdms_data_output);
    data_tdms_input = load(filename_tdms_data_input);
    % delta = data_tdms_input.delta;

    % num_measurement_cells_x = 2*ceil(0.5*length_measurement_cube/delta.x)+1;
    % num_measurement_cells_y = 2*ceil(0.5*length_measurement_cube/delta.y)+1;
    % num_measurement_cells_z = 2*ceil(0.5*length_measurement_cube/delta.z)+1;

    [num_measurement_cells_x, num_measurement_cells_y, num_measurement_cells_z] = size(data_tdms_input.ii);

    efield_tdms_output_raw = zeros(num_measurement_cells_x, num_measurement_cells_y, num_measurement_cells_z, 3);

    efield_tdms_output_raw(:) = data_tdms_output.campssample(:);
    efield_tdms_output = efield_tdms_output_raw(:, :, :, 1);
end
