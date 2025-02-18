dat = load('S014_converted_HT3D_0.mat');

%We are going to image with dat.data(:,b,c) as our z direction, dat.data(a,:,c) as our x direction, and dat.data(a,b,:) as our y direction.
xvec = dat.yvec;
yvec = dat.zvec;
zvec = dat.xvec;
pixel_size = dat.pixel_size;

%RI = zeros(numel(xvec),numel(yvec),numel(zvec));
%tic;
%for i=1:numel(xvec)
%    for j=1:numel(yvec)
%	for k=1:numel(zvec)
	    
%	    RI(i,j,k) = dat.data(k,i,j);
%	    if toc>1
%		fprintf(1,'(%d,%d,%d)\n',i,j,k);
%		tic;
%	    end
%	end
%   end
%end

RI = permute(dat.data,[2 3 1]);
RI = pagetranspose(RI);
outfile = 'S014_converted_HT3D_0_reshaped';
save(outfile,'RI','pixel_size','xvec','yvec','zvec','-V7.3');
