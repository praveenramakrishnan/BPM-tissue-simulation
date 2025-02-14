filename_tomocube_data_directory = '~/POL_NAS/tomocube_data/';
addpath(filename_tomocube_data_directory);
filename_tomocube_data = append(filename_tomocube_data_directory, 'S014_converted_HT3D_0_reshaped');
save(append(mfilename, '.mat'));
