This example simulates the propagation of light through the sphere based
model of the tissue. The randomly distributed ensemble of spheres of a specified concentration is generated in the volume.

# Input files
 - 'input_parameters_bpm.m' is used to specify the parameters related to the BPM simulation. 
 - 'input_parameters_sphere_refractive_index.m' is used to specify the parameters related to the refractive indices of the model.

# Running the simulation and generating the output
- Running 'main_bpm_sphere_ensemble.m' generates the outputs. The names of the generated output folders and files can be changed in the main file.
- The main file also contains a few simulation parameters.
- An output data file (by default 'output_data_sphere_code/') will store the data related to the centers and radii of the sphere ensemble.
- The output electric field data will be store in the specified output folder (by default starts with 'simulation_bpm_.../')
