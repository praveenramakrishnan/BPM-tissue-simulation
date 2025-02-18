function [wx,wy] = gauss_pol(th,ph, lambda, refind, NA_obj, eta);
    arguments
        th;
        ph;
        % Default parameters
        lambda = 920e-9;
        refind = 1.3333;
        NA_obj = 1.05;
        eta = 0.1; %relative magnitude of field at aperture
    end

    k = 2*pi/lambda;
    
    F = sqrt(log(1/eta));
      
    wx = exp( -(sin(th)*F/(NA_obj/refind)).^2 );
    wy = zeros(size(wx));

    % save gpvars;
end
