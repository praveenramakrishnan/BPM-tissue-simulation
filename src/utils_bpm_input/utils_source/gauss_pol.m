function [wx,wy] = gauss_pol(th,ph, lambda);
    if nargin < 3
        % Default value of wavelength
        lambda = 920e-9;
    end

    refind = 1.3333;
    NA_obj = 1.05;
    k = 2*pi/lambda;
    
    eta = 0.1; %relative magnitude of field at aperture
    F = sqrt(log(1/eta));
      
    wx = exp( -(sin(th)*F/(NA_obj/refind)).^2 );
    wy = zeros(size(wx));

    % save gpvars;
