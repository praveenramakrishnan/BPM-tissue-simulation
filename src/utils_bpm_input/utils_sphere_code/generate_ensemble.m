%[x,R] = generate_ensemble(x_min,x_max,z_min,z_max,rad_min,rad_max,N)
%
%Calculate an ensemble of spheres in the xz plane. All spheres
%centers must be contained within the region [x_min,x_max] x [y_min,y_max] x [z_min,zmax].
%
%The spheres have radii in the range [rad_min,rad_max]
%
%There will be N spheres
%
%Outputs:
%
%x - 3xN matrix of sphere centres
%R - 1xN matrix of radii
%
function [x,R] = generate_ensemble(x_min,x_max,y_min,y_max,z_min,z_max,rad_min,rad_max,N)

show_plot = false;

x=[];
R=[];

tic;
t=[];

while size(x,2)<N
    x_t=[rand(1)*(x_max-x_min-2*rad_max)+x_min+rad_max;
     rand(1)*(y_max-y_min-2*rad_max)+y_min+rad_max;   
	 rand(1)*(z_max-z_min-2*rad_max)+z_min+rad_max];
    
    R_t=rand(1)*(rad_max-rad_min)+rad_min;
    if size(x,2)>1
        if ~do_overlap_ensemble(R_t,x_t,R,x)
            x=[x x_t];
            R=[R R_t];
            if toc>5 % print progress every 5 seconds
                display('Generating sphere number: ');
                fprintf(1,'%d\n',size(x,2));
                tic;
            end
            %t=[t toc];
        end
    else
        x=[x x_t];
        R=[R R_t];
    end
    
end

if show_plot
    theta=linspace(0,2*pi,10);
    phi=linspace(-pi/2,pi/2,10);
    [phi,theta] = meshgrid(phi,theta);
    
    figure(1);clf;
    for i=1:numel(R)
        X=x(1,i)+R(i)*cos(phi).*cos(theta);
        Y=x(2,i)+R(i)*cos(phi).*sin(theta);
        Z=x(3,i)+R(i)*sin(phi);
        mesh(X,Y,Z);
        hold on;
    end
    axis equal;
end
