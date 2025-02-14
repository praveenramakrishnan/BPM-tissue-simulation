function [efield_propagated] = run_beam_propagation(efield_initial, num_grid_points_x, num_grid_points_y, ...
        num_grid_points_z, step_size_z, frequency_z_grid_nd, apply_phase_correction, num_grid_points_z_ri, ...
        refractive_index_data, lambda, refractive_index_average)
    
    efield_propagated = efield_initial;

    for iloop=1:num_grid_points_z 
        efield_propagated_homogeneous = bpm_propagate_homogeneous(efield_propagated, step_size_z, frequency_z_grid_nd);
        if apply_phase_correction
            iloop_refractive_index_data = ceil(iloop*num_grid_points_z_ri/num_grid_points_z);
            refractive_index_data_2d_mapped = map_refractive_index_grid( ...
                refractive_index_data(:,:, iloop_refractive_index_data), num_grid_points_x, num_grid_points_y);
            refractive_index_difference = refractive_index_data_2d_mapped - refractive_index_average;
            efield_propagated = bpm_apply_phase_correction(efield_propagated_homogeneous, step_size_z, refractive_index_difference, lambda);
        else 
            efield_propagated = efield_propagated_homogeneous;
        end
        display(append('Running iteration ', num2str(iloop), ', maximum field level = ', ...
            num2str(max(abs(efield_propagated(:))))));
    end

end

function [efield_propagated_homogeneous] = bpm_propagate_homogeneous(efield_initial, step_size_z, frequency_z_grid_nd)
    fft_efield_initial = fftshift(fft2(efield_initial));
    propagation_factor_homogeneous = exp(1j*2*pi*frequency_z_grid_nd*step_size_z);
    fft_efield_propagated_homogeneous = fft_efield_initial.*propagation_factor_homogeneous;
    efield_propagated_homogeneous = ifft2(ifftshift(fft_efield_propagated_homogeneous)); 
end

function [efield_propagated] = bpm_apply_phase_correction(efield_propagated_homogeneous, step_size_z, refractive_index_difference, lambda)
    phase_correction_factor = exp(1j*2*pi/lambda*refractive_index_difference*step_size_z);
    efield_propagated = efield_propagated_homogeneous.*phase_correction_factor;
end

function [refractive_index_data_2d_mapped] = map_refractive_index_grid(refractive_index_data_z, num_grid_points_x, num_grid_points_y)
    [num_points_data_x, num_points_data_y] = size(refractive_index_data_z);
    for iloop=1:num_grid_points_x
        for jloop=1:num_grid_points_y
            iloop_mapped = ceil(iloop*num_points_data_x/num_grid_points_x);
            jloop_mapped = ceil(jloop*num_points_data_y/num_grid_points_y);
            refractive_index_data_2d_mapped(iloop, jloop) = refractive_index_data_z(iloop_mapped, jloop_mapped);
        end
    end
end

