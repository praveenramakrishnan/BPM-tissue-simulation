function [efield_initial, x_grid, y_grid] = read_efield_input_data(filename_efield_input)
    data_efield_input = load(filename_efield_input);
    efield_initial = data_efield_input.efield_illumination{1}(:,:,1);
    x_grid = data_efield_input.x_grid;
    y_grid = data_efield_input.y_grid;
end

