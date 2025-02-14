function [wx,wy] = gauss_pol(th,ph);

    refind = 1.3333;
    lambda = 1300e-9;
    NA_obj = 1.05;
    k = 2*pi/lambda;
    
    eta = 0.1; %relative magnitude of field at aperture
    F = sqrt(log(1/eta));
      
    wx = exp( -(sin(th)*F/(NA_obj/refind)).^2 );
    wy = zeros(size(wx));

    save gpvars;
