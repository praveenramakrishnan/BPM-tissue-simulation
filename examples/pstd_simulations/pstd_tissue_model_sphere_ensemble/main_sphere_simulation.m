% Add path to library functions
clear all; close all;
run('../tdms_location.m');
addpath(append(tdms_root, 'tdms/tests/system/data/input_generation/matlab'));
addpath(genpath('../../../src/utils_simulation_setup/'));
addpath('../utils_sphere_tdms/');

% parameters
recompute_gridfile_background = false;
recompute_illumination_background = false;
recompute_tdms_setup_background = false;
recompute_background_simulation = false;

recompute_gridfile_sphere = true;
recompute_illumination_sphere = true;
recompute_tdms_setup_sphere = true;
recompute_sphere_simulation = true;

scattering_geometry = 0;
output_directory = append('artefacts_tdms_sphere_scattgeo_', num2str(scattering_geometry), '/');
if ~exist(output_directory)
    mkdir(output_directory);
end

display("Simulating scattering geometry " + scattering_geometry ...
    + " with output directory " + output_directory);
% Read input file
filename_input_pstd = 'input_file_sphere.m';

filename_input_variables_temp = append('variables_', strrep(filename_input_pstd, '.m', '.mat'));
filename_input_variables = append(output_directory, 'variables_input_file.mat');

filename_sphere_input = '../sphere_code/sphere_ensemble.mat';

% Create grid
[x_grid, y_grid, z_grid, lambda] = fdtd_bounds(filename_input_pstd);

if numel(y_grid) < 2
    y_grid = 0;
end

save(append(output_directory, 'variables_grid'), ...
    'x_grid', 'y_grid', 'z_grid','-v7.3');

movefile(filename_input_variables_temp, filename_input_variables);

% Parameters defining ensemble of spheres
input_params = load(filename_input_variables);
length_along_z = z_grid(end) - z_grid(1);
refractive_index_background = input_params.n_background;

if scattering_geometry==0
    radius_sphere_list = [length_along_z/10];
    x_center_sphere_list = [[0, 0, 0]];
    refractive_index_sphere_list = [refractive_index_background];
elseif scattering_geometry==1
    radius_sphere_list = [length_along_z/10];
    x_center_sphere_list = [[0, 0, 0]];
    refractive_index_sphere_list = [1.4];
elseif scattering_geometry==2
    radius_sphere_list = [length_along_z/10, length_along_z/10];
    x_center_sphere_list = [[0, 0, -length_along_z/4]; [0, 0, 0]];
    refractive_index_sphere_list = [1.4, 1.4];
elseif scattering_geometry==3
    radius_sphere_list = [length_along_z/10, length_along_z/10, length_along_z/10];
    x_center_sphere_list = [[0, 0, -length_along_z/4]; ...
        [0, 0, length_along_z/4];
        [0, 0, 0]];
    refractive_index_sphere_list = [1.4, 1.4, 1.4];
elseif scattering_geometry==4
    input_variables_sphere = load(filename_sphere_input);
    radius_sphere_list  = [input_variables_sphere.radius_sphere];
    x_center_sphere_list = [input_variables_sphere.x_sphere_center]';
    refractive_index_sphere_list = 1.39*ones(size(radius_sphere_list));
end

save(append(output_directory, 'variables_sphere.mat'), 'radius_sphere_list', ...
    'x_center_sphere_list', 'refractive_index_sphere_list', 'refractive_index_background', '-v7.3');

% 1) Set up executable file for background material simulation

% Background material properties
filename_gridfile_background = append(output_directory, 'gridfile_background.mat');
filename_illumination_background = append(output_directory, 'illumination_field_background.mat');
filename_setup_background = append(output_directory,  'setup_field_background.mat');
if recompute_gridfile_background
    display('Begin creating gridfile for background');
    composition_matrix = [];
    material_matrix = [];
    save(filename_gridfile_background, 'material_matrix', 'composition_matrix', '-v7.3');
    display('End creating gridfile for background');
end

if recompute_illumination_background
    display('Begin iterate fdtd calculation for background illumination');
    % use iteratefdtd_matrix to set up illumination file
    iteratefdtd_matrix(filename_input_pstd, 'illsetup', ...
        filename_illumination_background, filename_gridfile_background, '');
    display('End iterate fdtd calculation for background illumination');
end

if recompute_tdms_setup_background
    display('Begin iterate fdtd calculation for background file setup');
    % use iteratefdtd_matrix to set up file for tdms execution
    iteratefdtd_matrix(filename_input_pstd, 'filesetup', filename_setup_background,...
        filename_gridfile_background, filename_illumination_background);
    display('End iterate fdtd calculation for background file setup');
end

% % 2) Set up executable file for sphere simulation 
filename_gridfile_sphere = append(output_directory, 'gridfile_sphere.mat');
filename_illumination_sphere = append(output_directory, 'illumination_field_sphere.mat');
filename_setup_sphere = append(output_directory, 'setup_field_sphere.mat');
filename_refractive_index_data_3d = append(output_directory, 'refractive_index_data_3d.mat');
filename_refractive_index_data_2d = append(output_directory, 'refractive_index_data_2d.mat');

if recompute_gridfile_sphere
    display('Begin creating gridfile for spheres');
    [composition_matrix, material_matrix] = create_material_grid(filename_input_pstd, ...
        radius_sphere_list, x_center_sphere_list, refractive_index_sphere_list);

    refractive_index_ndgrid = create_refractive_index_sphere(filename_input_pstd, ...
        radius_sphere_list, x_center_sphere_list, ...
        refractive_index_sphere_list, refractive_index_background);

    save(filename_gridfile_sphere, 'material_matrix', 'composition_matrix', '-v7.3');
    save(filename_refractive_index_data_3d, 'refractive_index_ndgrid', '-v7.3');

    display('End creating gridfile for spheres');
end

if recompute_illumination_sphere
    display('Begin iterate fdtd calculation for sphere illumination');
    % use iteratefdtd_matrix to set up illumination file
    iteratefdtd_matrix(filename_input_pstd, 'illsetup', ...
    filename_illumination_sphere, filename_gridfile_sphere, '');
    display('End iterate fdtd calculation for sphere illumination');
end

if recompute_tdms_setup_sphere
    display('Begin iterate fdtd calculation for sphere file setup');
    % use iteratefdtd_matrix to set up file for tdms execution
    iteratefdtd_matrix(filename_input_pstd, 'filesetup', filename_setup_sphere,...
        filename_gridfile_sphere, filename_illumination_sphere);
    display('End iterate fdtd calculation for sphere file setup');
end

% Run the tdms
filename_output_background =  append(output_directory, 'output_fields_background.mat');
filename_output_sphere =  append(output_directory, 'output_fields_sphere.mat');
filename_tdms_iterations_background = append(output_directory, 'output_tdms_iterations_background.txt');
filename_tdms_iterations_sphere = append(output_directory, 'output_tdms_iterations_sphere.txt');

if recompute_background_simulation
    display('Begin tdms simulation for background');
    run_tdms_command_background = append('tdms ', filename_setup_background, ' ', ...
        filename_output_background, ' >> ', ...
        filename_tdms_iterations_background);
    system(run_tdms_command_background);
    display('End tdms simulation for background');
end

if recompute_sphere_simulation
    display('Begin tdms simulation for sphere');
    run_tdms_command_sphere = append('tdms ', filename_setup_sphere, ' ', ...
        filename_output_sphere, ' >> ', ...
        filename_tdms_iterations_sphere);
    system(run_tdms_command_sphere);
    display('End tdms simulation for sphere');
end

% Clean up directory and move all outputs to right folder
if ~isempty(dir('*.mat'))
    system(append('mv *.mat ', output_directory));
end
