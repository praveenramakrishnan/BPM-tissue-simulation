%function [E] = efield_focused_gauss(X,Y,Z,delta);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%
function [E] = efield_focused_gauss(X,Y,Z, lambda)
    if nargin < 4
        % Default value of wavelength
        lambda = 920e-9;
    end
 
    E{1} = zeros(size(X));
    E{2} = zeros(size(X));
    E{3} = zeros(size(X));

    vertices = [X(:) Y(:) Z(:)];
    nvec = 1.3333; 
    hvec = [];
    NA = 1.05;
    ntheta = 1600;
    nphi = 2000;
    %ntheta = 200;
    %nphi = 200;
    
    polfun = @gauss_pol;
    
    [Ep,Em] = focstratfield_general_pol(vertices,nvec,hvec,NA,lambda,ntheta,nphi,polfun);
    [Ep0,Em0] = focstratfield_general_pol([0 0 0],nvec,hvec,NA,lambda,ntheta,nphi,polfun);
    E{1}(:) = Ep(:,1)/Ep0(1,1);
    E{2}(:) = Ep(:,2)/Ep0(1,1);

