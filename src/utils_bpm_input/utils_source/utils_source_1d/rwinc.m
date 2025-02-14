%function [E,d] = rwinc(NA,k,xlims,ylims,zlims,nintegral)
%
%NA - numerical aperture
%k - wave number
%xlims - [xmin xmax Nx]
%ylims - [ymin ymax Ny]
%zlims - [zmin zmax Nz]
%nintegral - number of points in abscissa
%
%Returns - 
%
%E - A cell array representing the electric field on the image plane arranged as:
%       E{1} - E field in x direction
%       E{2} - E field in y direction
%       E{3} - E field in z direction
%
%d - A cell array representing the optical coordinates on the image plane arranged as:
%       d{1} - x optical coordinates
%       d{2} - y optical coordinates
%       d{3} - z optical coordinates
%
%For example:
%
%[E,d] = rwinc(.85,2*pi/405e-9,[-.5e-6 .5e-6 101],[-.5e-6 .5e-6 101],[0 0 1],50);
%
%I = E{1}.*conj(E{1})+E{2}.*conj(E{2})+E{3}.*conj(E{3});
%h = imagesc(d{1},d{2},real(transpose(I)));
%shading interp;
%axis xy;
%axis equal;
%axis tight;
%colormap(gray);

function [E,d] = rwinc(NA,k,xlims,ylims,zlims,nintegral)
    
    x = linspace(xlims(1),xlims(2),xlims(3));
    y = linspace(ylims(1),ylims(2),ylims(3));
    z = linspace(zlims(1),zlims(2),zlims(3));
    
    d = {x,y,z};
    
    [X,Y,Z] = ndgrid(x,y,z);
    
    [m,n,o] = size(X);
    
    xvec = reshape(X,1,m*n*o);
    yvec = reshape(Y,1,m*n*o);
    zvec = reshape(Z,1,m*n*o);
        
% $$$     NA - numerical aperture of focusing lens
% $$$     k  - wavenumber
% $$$     xverts - x component of vertices
% $$$     yverts - y component of vertices
% $$$     zverts - z component of vertices
% $$$     nintegral - number of integration steps
    
    [Ex,Ey,Ez] = richwolfincident(NA,k,xvec,yvec,zvec,nintegral);
    
    E = {reshape(Ex,m,n,o),reshape(Ey,m,n,o),reshape(Ez,m,n,o)};
    
    
