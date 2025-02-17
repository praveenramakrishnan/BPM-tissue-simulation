filename_tomocube_data_directory = '~/POL_NAS/tomocube_data/';
addpath(filename_tomocube_data_directory);
filename_tomocube_data = append(filename_tomocube_data_directory, 'S014_converted_HT3D_0_reshaped');
% If not known a priori, the background refractive index may be computed as the average of
% the refractive index data
refractive_index_background = 1.3333; 
save(append(mfilename, '.mat'));
