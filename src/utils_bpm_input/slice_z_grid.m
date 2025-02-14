function [z_grid_slice] = slice_z_grid(z_grid_full, num_slices, slice_index)
    % Read input parameters
    K = length(z_grid_full);
    K_start = max(1, (slice_index-1)*round(K/num_slices));
    K_end = slice_index*round(K/num_slices);

    % Define grid
    z_grid_slice = z_grid_full(K_start:K_end);
end
