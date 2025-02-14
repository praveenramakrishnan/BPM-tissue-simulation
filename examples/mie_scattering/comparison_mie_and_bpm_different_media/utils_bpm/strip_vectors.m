function [vector_stripped] = strip_vectors(input_vector, Nx, Ny)
    [Nx_full, Ny_full] = size(input_vector);
    vector_stripped = input_vector(ceil((Nx_full-Nx)/2)+1:floor((Nx_full+Nx)/2), ceil((Ny_full-Ny)/2)+1:floor((Ny_full+Ny)/2), :);
end
