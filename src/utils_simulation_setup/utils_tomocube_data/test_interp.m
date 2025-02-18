dat=load('S014_converted_HT3D_0_reshaped');

xi = linspace(100e-6,250e-6,500);
%yi = linspace(0,50e-6,100);
yi = 67.233e-6;
zi = linspace(50e-6,200e-6,500);

[XI,YI,ZI] = ndgrid(xi,yi,zi);

%RI is in the meshgrid format, which is why it has size 154x3860x3860
[RIint] = interp_data(dat.RI,dat.xvec,dat.yvec,dat.zvec,[XI(:) YI(:) ZI(:)]);
RI = zeros(numel(xi),numel(yi),numel(zi));
RI(:) = RIint;
RI2 = pagetranspose(dat.RI);


figure(1);clf;
imagesc(dat.zvec,dat.xvec,squeeze(RI2(:,72,:)));colormap gray;axis equal;
ca=caxis;

figure(2);clf;
imagesc(zi,xi,squeeze(RI));colormap gray;axis equal;
ax=axis;
caxis(ca);

figure(1);axis(ax);

%now test varying yi
yi =  linspace(0,100e-6,250);
[XI,YI,ZI] = ndgrid(xi,yi,zi);

T = dat.yvec(90)-dat.yvec(54);%useful range
y0 = dat.yvec(54);
[RIint] = interp_data(dat.RI,dat.xvec,dat.yvec,dat.zvec,[XI(:) sawtooth(YI(:)-y0,T)+y0 ZI(:)]);
RI = zeros(numel(xi),numel(yi),numel(zi));
RI(:) = RIint;



figure(3);clf;
imagesc(yi,xi,squeeze(RI(:,:,250)));colormap gray;axis equal;
