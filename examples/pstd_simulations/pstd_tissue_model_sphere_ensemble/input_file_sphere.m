% The inputs are divided into the following sections:
% 1. Grid,
% 2. Source,
% 3. Simulation type,
% 4. FDTD specific,
% 5. Output,
% 6. Perfectly Matched Layer (PML),

% 0) Basic parameters and constants

c_speed_of_light = 2.997924580105029e+08;
lambda = 1300e-9;

%Use linear interpolation
use_bli=0;

%Select FDTD or PSTD simulation
use_pstd=1;

% 1) Grid
% These control the resolution and size of the grid as well as the
% default material in which scattering objects are embedded. Parameters
% to set the multilayer structure are included in this category.

% This specifies specifies the the dimensionality of the simulation.
% It can can take the values 'TE’, 'TM’ and '3’ which represent the TE, TM and
% three dimensional simulations simulations respectively. The default value is 3 dimension. 
dimension = '3';

length_along_x = 150*lambda;
length_along_y = 150*lambda;
length_along_z = 55*lambda;

% size of Yee cell in metres
% set grid size to be same as pixel size of refractive index data
delta.x = lambda/10; 
delta.y = lambda/10;
delta.z = lambda/10;

% Define the grid size
% I,J, K:The number of non-PML cells in the FDTD grid in in the the x, y, z directions.
% I = 640;
% J = 640;
% K = 170;

I = round(length_along_x/delta.x);
J = round(length_along_y/delta.y);
K = round(length_along_z/delta.z);

% Each element of multilayer gives the index of the Yee cell along the z-direction of an interface.
% A simulation without a multilayer structure, ie, homogeneous space would have multilayer set to [].
multilayer = [];

% Possibly complex non-dispersive relative permittivity of the background
% material which may be a multilayer structure.
n_background = 1.3333;
epsr = [n_background^2]; %epsr = (refractive index)^2

% radian plasma frequency of each layer
%wp_vec = [0];

% collision frequency of each layer
%vc_vec = [0];

% relative permeability of the background material. 
mur = 1;

% do not interpolate material properties
intmatprops = 0;

% piecewise linear profile of the multilayer structure when looking at the xz plane.
% structure = [];


% This is the name of the file which describes the scattering object
% within the FDTD grid.
% material_file = '';

% 2) Source
% These determine the properties of the incident electromagnetic field.

% frequency in Hz
f_an =  c_speed_of_light/lambda; 

% describe when the incident field is introduced into the grid
source_interface_location_z = 10;
interface.I0 = [5 0];
interface.I1 = [I-5 0];
interface.J0 = [5 0];
interface.J1 = [J-5 0];
interface.K0 = [source_interface_location_z 1];
interface.K1 = [K-5 0];

% these are the function names used to generate the field
efname = 'efield_plane';
hfname = '';

% this defines the point about which the illumination is centred in
% the so called 'interior' coordinate system
% this means the the illumination is actually focussed on a point 25
% cells in front of the interface.
% K_focus = K-50;
length_measurement_cube = 10e-6;
K_focus = K-ceil(0.5*length_measurement_cube/delta.z);
illorigin = [ceil(I/2), ceil(J/2), K_focus];
% +- 2 microns around origin (4 micrometer cube).

% the z coordinate of illorigin cell
% z_launch = 0;

%the wavelength width (in m). This corresponds to the FWHM of the
wavelengthwidth = 100e-9;

%Not a compact source condition
% compactsource=1;
compactsource=1;

%this is the width of the pulse in terms of the time step
%gpulsewidth = 400;


% k_vec = 2*pi/lambda;

% f_ex_vec=2.997924580105029e+08*k_vec/2/pi;

exdetintegral=0;
k_det_obs=5;
NA_det=(7e-3/2)/36e-3;
beta_det=25/36;
detmodevec=1;
detsensefun='gaussian_d_telesto';

% 3) Simulation type
% Determines what kind of simulation is performed.
% Allows, for example, a simulation type amenable to analysis of the
% field in the time domain or the frequency domain to be chosen.

%this is the kind of source mode, can be 'pulsed' or 'steadystate'. 
sourcemode = 'pulsed';

%this defines the run mode of the simulation, can be 'analyse' or
%'complete'. 'analyse' means that sub results can be saved using
%the statements in outputs_array. When complete is specified, only
%the final results will be saved using the outputs_array statements.
%runmode = 'analyse';
runmode = 'complete';

% 4) FDTD specific
% Parameters which are specific to the FDTD algorithm.
% These must be chosen by someone with some experience with FDTD.
% Some of the source parameters also fit into this category.


%time step - subject to restriction
factor_time_step_size = 1;
dt = 1/factor_time_step_size*2/sqrt(3)/pi*delta.x/(3e8/1.3)*.95;

%define the number of time steps
Nt=4000*factor_time_step_size;

k_vec = 2*pi/lambda;

f_ex_vec=asin(k_vec*c_speed_of_light*dt/2)/(pi*dt);

% 5) Output
% Determines which field parameters are calculated and outputed by the algorithm.

%  enables various field components at selected points in the grid to be output at
% each iteration of the simulation. In general this should be set as outputs_array={}
% as this option may only be used when runmode is set to 'analyse’
% not used as run mode is complete
outputs_array ={};    

% calculate phasors in the entire volume of the FDTD grid if set to 1.
exphasorsvolume = 0;

%this determines whether or not to extract phasors around a
%specified surface
exphasorssurface = 0;

%this specifies a surface to extract the phasors at. These
%quantities are in interior coordinate system;
%has the form [I0 I1 J0 J1 K0 K1] which defines the extremes of a
%cuboid wihch defines the surface to extract phasors at
%These should be set so that the interpolation scheme can work
phasorsurface = [5 I-5 5 J-5 5 K-5];

% divide FDTD lattice dimension
% phasorinc = [1, 1, 1,];

dI = ceil(0.5*length_measurement_cube/delta.x);
dJ = ceil(0.5*length_measurement_cube/delta.y);
dK = ceil(0.5*length_measurement_cube/delta.z);

% [ii,jj,kk] = ndgrid(ceil(I/2),ceil(J/2),1:K);
% K-20-dk: k-20+dk  dk = ceil(4 microns/delta.z) 
[ii,jj,kk] = ndgrid(ceil(I/2)-dI:ceil(I/2)+dI,ceil(J/2)-dJ:ceil(J/2)+dJ,K_focus-dK:K_focus+dK);

% [ii,jj,kk] = ndgrid(3,3,1:K);
% [ii,jj,kk] = ndgrid(300:(I-300),300:(J-300),[50:50:300 (illorigin(3)-50):(illorigin(3)+50)]);
campssample.vertices = [ii(:) jj(:) kk(:)];
campssample.components = [1 2 3];

% fieldsample.i = 300:(I-300);
% fieldsample.j = 300:(J-300);
% fieldsample.k = [50:50:300 (illorigin(3)-50):(illorigin(3)+50)];
% fieldsample.n = [2 4];

fieldsample.i = [];
fieldsample.j = []; 
fieldsample.k = []; 
fieldsample.n = []; 

% 6) Perfectly matched later
% This is used to model unbounded free space scattering with a finite grid.
% These must be set by an experienced user.

%order of the PML conductivity profile curve
n = 4;

%maximum reflection at PML
R0 = 1e-7;

% kappa parameter in PML useful to terminate conductive and dispersive materials
% set to 1 for terminating dielectrics 
kappa_max = [1];

%number of PML cells in each direction, this layer absorbs light
%but does not reflect it
Dxl = 10;
Dxu = 10;
Dyl = 10;
Dyu = 10;
Dzl = 10;
Dzu = 10;

save(append('variables_', mfilename, '.mat'));
