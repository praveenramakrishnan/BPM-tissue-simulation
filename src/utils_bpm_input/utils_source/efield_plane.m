%function [E] = efield_focused_gauss(X,Y,Z,delta);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%
function [E] = efield_plane(X,Y,Z, lamdba, nvec)
    if nargin < 5
        % Default value of wavelength
        lambda = 920e-9;
        nvec = 1.3333;
    end
    lambda_medium = lambda/nvec; % wavelength in the medium

    E{1} = exp(sqrt(-1)*2*pi/lambda_medium*Z);
    E{2} = zeros(size(X));
    E{3} = zeros(size(X));
