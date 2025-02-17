This codebase contains the implementation of the beam propagation method (BPM) for
propagation  of light through brain tissue models.
Different examples are provided to demonstrate the usage.
The method can be used to simulate large tissue models.

# Solver methods
## 1. Beam Propagation Method (BPM)
A fast fourier transform (FFT) based BPM is implemented. The computation is carried out
with two planes at a time. The electric field incident on a plane is give as the input
which is split into its constituent plane wave components and propagated to the output
plane and recombined. The various parameters of the solver are provided in the examples
folder. 

## 2. Mie solver
The code for Mie solver is provided which can be used to validate the BPM solver for a
single spherical scatterer. Mie theory provides a series solution for the electromagnetic
field due to plane wave illumination on a sphere. The code provided assumes that the
sphere is centered at the origin. The incident wave is assumed to propagate along the
z-axis and linearly polarized along the x-axis. 

## 3. Pseudo Spectral Time Domain (PSTD) Method
The PSTD is a numerical method to solve Maxwell's equation. While this method can be
considered as more rigous than BPM, it is computationally demanding. We provide examples
where we compare the solution obtain using both PSTD and BPM for small enough problems.
The PSTD method is implemented in the [TDMS repository](https://github.com/UCL/TDMS) and
we provide examples of the its usage for the problem of interest.

# Source Field Computation
The BPM can be run using a variety of illumination sources including plane wave and
focused beam. The code for these two sources are provided.

The focused beam illumination is computed for a partially filled aperture. The method for
computing the focused illumination is described in the paper: ["Peter RT Munro. *Tool for
simulating the focusing of arbitrary vector beams in free-space and stratified media*.
In: Journal of biomedical optics 23.9 (2018), pp.
090801–090801"](https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-23/issue-09/090801/Tool-for-simulating-the-focusing-of-arbitrary-vector-beams-in/10.1117/1.JBO.23.9.090801.full) 

# Tissue Models
The refrctive index of the propagation needs to be provided to carry out the simulations.
The following two approaches are used in the provided examples.
## 1. Discrete model with ensemble of spheres
The tissue can be modelled using an ensemble of spheres randomly distributed in the
background medium. The size and concentration of spheres can be computed using the Mie
theory to match the required required scattering parameters such as the scattering
coefficient and average cosine of phase function. [Open source tools](https://omlc.org/calc/index.html) are
available for the calculation.
## 2. Measured refractive index data
One can use refractive index data for the tissue obtained using actual measurements.
The tomocube data set used in some of the examples which provides the refractive index of
mouse brain tissue. As mentioned below, this data set is not openly available yet.

# Installation
The main BPM is has been tested on a Linux machine with Matlab 2021.
The codebase is mostly in Matlab while a few functions are written in c/c++ and compiled
as mex files. Many examples do not have any other requirement.

A few examples use the PSTD method which require installing
the [TDMS repository](https://github.com/UCL/TDMS).
The link to the folder can be updated in the input file provided in the folders containig
the PSTD example.

Running the example with tomocube data requies access to the dataset described in the paper
["Lee, Ariel J., et al.  *Volumetric Refractive Index Measurement and Quantitative Density Analysis of Mouse Brain Tissue with Sub‐Micrometer Spatial Resolution*. Advanced Photonics Research 4.10 (2023): 230011".](https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/adpr.202300112)
This data set is currently not publicly available but might be requested from the authors of the paper. 
