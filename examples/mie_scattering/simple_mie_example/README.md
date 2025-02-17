This program computes the electric field due to a spherical scatterer centered at origin
under plane wave illumination.

The function "mie_series" calculates the complex amplitude of the electric and magnetic
fields for a plane wave, linearly polarised in the x-direction, propagating in the
positive z-direction, for a sphere centered at the origin.

Where:

- Nterms: The number of terms to take in the Mie series, can be set
         by when convergence occurs or using the Wiscombe criterion (https://opg.optica.org/ao/fulltext.cfm?uri=ao-19-9-1505&id=23949).
- incident: Include incident wave if set to 1
- n_sphere: The refractive index of the scattering sphere
- n_free: The refractive index of the material the sphere is embedded in
- E0: complex amplitude of incident Wave
- sphere_radius: Radius of the sphere
- lambda: Wavelength of incident illumination in a vacuum
- vertices: Matrix of vertices of points to evaluate field at, has
           dimension Nv x 3, where Nv is the number of vertices

- Eout: Matrix of electric field complex amplitudes at each
vertice, has same dimension as vertices.
- Hout: Matrix of magnetic field complex amplitudes at each
vertice, has same dimension as vertices.
