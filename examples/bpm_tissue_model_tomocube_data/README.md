This example simulates the propagation of light using the measured refractive index of
moise brain tissue (tomocube data).

# Input files
 - 'input_parameters_bpm.m' is used to specify the parameters related to the BPM simulation. 
 - 'input_parameters_tomocube_data.m' is used to specify the filename of the tomocube data. This data set is currently not publicly available but might be requested from the authors of the paper "Lee, Ariel J., et al.  *Volumetric Refractive Index Measurement and Quantitative Density Analysis of Mouse Brain Tissue with Sub‐Micrometer Spatial Resolution*. Advanced Photonics Research 4.10 (2023): 230011.". 

# Running the simulation and generating the output
- Running 'main_bpm_sphere_ensemble.m' generates the outputs. The names of the generated output folders and files can be changed in the main file.
- The main file also contains a few additional simulation parameters.
  For example, the variable source_type can be used to specify either
  a plane wave or focused beam source.
  Alternatively, the source can be specified directly by specifying
  the function 'efield_illumination_function'.
- An output data file (by default 'output_data_sphere_code/') will store the data related to the centers and radii of the sphere ensemble.
- The output electric field data will be store in the specified output folder (by default has the form 'simulation_bpm_.../')
