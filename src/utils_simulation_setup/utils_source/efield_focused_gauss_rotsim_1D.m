function [E] = efield_focused_gauss_rotsim_1D(X,Y,Z, lambda, nvec, NA, ntheta)
    arguments
        X; 
        Y;
        Z;
        % Default parameters
        lambda = 920e-9; % wavelength
        nvec = 1.3333; % background refractive index
        NA = 1.05; % numerical aperture
        ntheta = 10000; % number of integration points
    end

    RHO = sqrt(X.^2 + Y.^2);
    rho_vec = linspace(0,max(RHO(:)),ceil(max(RHO(:))/lambda*50));
 
    E{1} = zeros(size(X));
    E{2} = zeros(size(X));
    E{3} = zeros(size(X));

    vertices = [X(:) Y(:) Z(:)];
    hvec = [];
    
    theta = pi/4;
    x = rho_vec*cos(theta);
    y = rho_vec*sin(theta);
    z = ones(size(x))*Z(1);
    vertices = [x(:) y(:) z(:)];
    
    eta = 0.1; %relative magnitude of field at aperture
    F = sqrt(log(1/eta));

    [Ex,Ey,Ez] = richwolfincident(NA/nvec,2*pi/lambda*nvec,x,y,z,F,ntheta);
    
    I0 = Ex;%Ep(:,1);
    I2 = Ey;%Ep(:,2);
    
    THETA = atan2(Y,X);
    
    I0_int = zeros(size(X));
    %I1_int = zeros(size(X));
    I2_int = zeros(size(X));
    
    I0_int(:) = interp1(rho_vec,I0,RHO(:),'spline',0);
    %I1_int(:) = interp1(rho_vec,I1,RHO(:));
    I2_int(:) = interp1(rho_vec,I2,RHO(:),'spline',0);
    
    Ex = I0_int + cos(2*THETA).*I2_int;
    Ey = sin(2*THETA).*I2_int;
    
    [Ex0,Ey0,Ez0] = richwolfincident(NA/nvec,2*pi/lambda*nvec,0,0,0,F,ntheta);
    if abs(Ex0)>0
        Ex = Ex/Ex0;
    end

    if abs(Ey0)>0
        Ey = Ey/Ey0;
    end

    E{1}(:) = Ex;
    E{2}(:) = Ey;

    % save efgr_vars;
end
