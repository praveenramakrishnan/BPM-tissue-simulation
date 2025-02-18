%vertsi = Nx3
%
function [RIint] = interp_data(RI,xvec,yvec,zvec,vertsi);

    [X,Y,Z] = meshgrid(xvec,yvec,zvec);
    RIint = interp3(X,Y,Z,RI,vertsi(:,1),vertsi(:,2),vertsi(:,3),'cubic');
