infile = '230602.155856.mouse brain tissue.014.brain tissue.A1.S014_converted_HT3D_0.raw';
outfile = 'S014_converted_HT3D_0';
data = flipud(raw2mat('230602.155856.mouse brain tissue.014.brain tissue.A1.S014_converted_HT3D_0.raw', []));

pixel_size = [0.155433 0.155433 0.946946]*1e-6;
xvec = (0:(size(data,1)-1))*pixel_size(1);
yvec = (0:(size(data,2)-1))*pixel_size(2);
zvec = (0:(size(data,3)-1))*pixel_size(3);

save(outfile,'data','pixel_size','xvec','yvec','zvec','-V7.3');
