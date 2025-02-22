function [E] = efield_focused_lambda_1300nm(X, Y, Z)
    lambda = 1300e-9; % wavelength
    nvec = 1.3333; % background refractive index
    NA = 1.05; % numerical aperture
    ntheta = 2000; % number of integration points
    E = efield_focused_gauss_rotsim_1D(X, Y, Z, lambda, nvec, NA, ntheta);
end
