%function [E] = efield_focused_gauss(X,Y,Z,delta);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%
function [E] = efield_focused_gauss(X,Y,Z, lambda, nvec, NA, ntheta, nphi, polfun)
    arguments
        X;
        Y;
        Z;
        % Default parameters
        lambda = 920e-9;
        nvec = 1.3333; 
        NA = 1.05;
        ntheta = 1600;
        nphi = 2000;
        polfun = @gauss_pol;
    end
 
    E{1} = zeros(size(X));
    E{2} = zeros(size(X));
    E{3} = zeros(size(X));

    vertices = [X(:) Y(:) Z(:)];
    hvec = [];
    
    [Ep,Em] = focstratfield_general_pol(vertices,nvec,hvec,NA,lambda,ntheta,nphi,polfun);
    [Ep0,Em0] = focstratfield_general_pol([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,polfun);
    E{1}(:) = Ep(:,1)/Ep0(1,1);
    E{2}(:) = Ep(:,2)/Ep0(1,1);
