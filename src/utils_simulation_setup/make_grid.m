function [x_grid, y_grid, z_grid, lambda] = make_grid(filename_input_parameters)
    % Read input parameters
    data_input = load(filename_input_parameters);
    lambda = data_input.lambda;
    I = data_input.I;
    J = data_input.J;
    K = data_input.K;
    delta = data_input.delta;
    illorigin = data_input.illorigin;
    z_launch = data_input.z_launch;

    % Define grid
    x_grid = ((1:I) - illorigin(1))*delta.x;
    y_grid = ((1:J) - illorigin(2))*delta.y;
    z_grid = ((1:K) - illorigin(3))*delta.z + z_launch;
end
