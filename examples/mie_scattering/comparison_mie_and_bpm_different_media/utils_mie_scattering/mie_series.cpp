#include "math.h"
#include "mex.h"
#include <complex>
#include "matrix.h"

using namespace std;
#include "matlabio.h"
#include "elem_functions.h"

#define N_IN 8
#define N_OUT 2

/*
Expects:

Inputs:
-Nterms: The number of terms to take in the Mie series
-incident: Include incident wave if set to 1
-n_sphere: The refractive index of the scattering sphere
-n_free: The refractive index of the material the sphere is embedded in
-Eo: complex amplitude of incident wave
-sphere_radius: Radius of the sphere
-lambda: Wavelength of incident illumination in a vacuum
-vertices: Matrix of vertices of points to evaluate field at

Outputs:
E
H

*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  
  char t_buf[100];
  int Nterms, Npoints, count = 0;
  complex<double> n_sphere, n_free, Eo;
  double sphere_radius, lambda;
  complex<double> *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
  double **Er, **Ei, **Hr, **Hi, **vertices, *theta_vec, *rad_vec, *phi_vec;
  int ndim, incident;

   mwSize *dims;dims = (mwSize *)malloc(3*sizeof(mwSize));
   
   const mwSize *dimptr;

  if(nlhs != N_OUT ){
    sprintf(t_buf,"%d inputs required",N_IN);
    mexErrMsgTxt(t_buf);
  }

  if(nrhs != N_IN){
    sprintf(t_buf,"%d outputs required",N_OUT);
    mexErrMsgTxt(t_buf);
  }
    

  Nterms = (int)*mxGetPr(prhs[count++]);
  incident = (int)*mxGetPr(prhs[count++]);
  
  if(mxIsComplex(prhs[count]))
    n_sphere = complex<double>(*mxGetPr(prhs[count]),*mxGetPi(prhs[count]));
  else
    n_sphere = complex<double>(*mxGetPr(prhs[count]),0.);
  count++;

  if(mxIsComplex(prhs[count]))
    n_free = complex<double>(*mxGetPr(prhs[count]),*mxGetPi(prhs[count]));
  else
    n_free = complex<double>(*mxGetPr(prhs[count]),0.);
  count++;
  if(mxIsComplex(prhs[count]))
    Eo = complex<double>(*mxGetPr(prhs[count]),*mxGetPi(prhs[count]));
  else
    Eo = complex<double>(*mxGetPr(prhs[count]),0.);
  count++;
  sphere_radius = *mxGetPr(prhs[count++]);
  lambda = *mxGetPr(prhs[count++]);
  ndim = mxGetNumberOfDimensions(prhs[count]);
  if(ndim != 2)
    mexErrMsgTxt("The vertices matrix should be two dimensional");
  dimptr =  mxGetDimensions(prhs[count]);
  Npoints = dimptr[0];
  if(dimptr[1]!=3){
    fprintf(stderr,"%d %d\n",dimptr[0],dimptr[1]);
    mexErrMsgTxt("The vertices matrix should have dimension Npoint x 3");
  }
  
  vertices = castMatlab2DArray(mxGetPr(prhs[count]), Npoints, 3);
  
  /*
  fprintf(stderr,"Nterms: %d\n",Nterms);
  fprintf(stderr,"n_sphere: %e+i%e\n",real(n_sphere),imag(n_sphere));
  fprintf(stderr,"n_free: %e+i%e\n",real(n_free),imag(n_free));
  fprintf(stderr,"Eo: %e + i%e\n",real(Eo), imag(Eo));
  fprintf(stderr,"sphere: radius: %e\n",sphere_radius);
  fprintf(stderr,"lambda: %e\n",lambda);
  fprintf(stderr,"N vertices: %d\n",Npoints);
  */

  ndim = 2;
  dims[0] = Npoints;
  dims[1] = 3;

  plhs[0] = mxCreateNumericArray( ndim, (const mwSize*) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  plhs[1] = mxCreateNumericArray( ndim, (const mwSize*) dims, mxDOUBLE_CLASS, mxCOMPLEX);

  rad_vec   = (double *)malloc(Npoints*sizeof(double));
  theta_vec = (double *)malloc(Npoints*sizeof(double));
  phi_vec   = (double *)malloc(Npoints*sizeof(double));
    
  Er = castMatlab2DArray(mxGetPr(plhs[0]), Npoints, 3);
  Ei = castMatlab2DArray(mxGetPi(plhs[0]), Npoints, 3);
  Hr = castMatlab2DArray(mxGetPr(plhs[1]), Npoints, 3);
  Hi = castMatlab2DArray(mxGetPi(plhs[1]), Npoints, 3);

  Ex = (complex<double> *)malloc(Npoints*sizeof(complex<double>));
  Ey = (complex<double> *)malloc(Npoints*sizeof(complex<double>));
  Ez = (complex<double> *)malloc(Npoints*sizeof(complex<double>));

  Hx = (complex<double> *)malloc(Npoints*sizeof(complex<double>));
  Hy = (complex<double> *)malloc(Npoints*sizeof(complex<double>));
  Hz = (complex<double> *)malloc(Npoints*sizeof(complex<double>));

  /*First generate the polar coordinates*/
  double x, y, z;
  for(int i=0;i<Npoints;i++){
    x = vertices[0][i];
    y = vertices[1][i];
    z = vertices[2][i];
    
    rad_vec[i]   = sqrt(x*x + y*y + z*z);
    theta_vec[i] = atan2(sqrt(x*x + y*y),z);
    phi_vec[i]   = atan2(y,x);

  }

  //calculate the field
  eval_elem_functions(Nterms, incident, n_sphere, n_free, sphere_radius, lambda, 
		      rad_vec, theta_vec, phi_vec, Npoints, Eo,
		      Ex, Ey, Ez, 
		      Hx, Hy, Hz);

  //populate the field vectors
  for(int i=0;i<Npoints;i++){
    Er[0][i] = real(Ex[i]);
    Ei[0][i] = imag(Ex[i]);
    Er[1][i] = real(Ey[i]);
    Ei[1][i] = imag(Ey[i]);
    Er[2][i] = real(Ez[i]);
    Ei[2][i] = imag(Ez[i]);

    Hr[0][i] = real(Hx[i]);
    Hi[0][i] = imag(Hx[i]);
    Hr[1][i] = real(Hy[i]);
    Hi[1][i] = imag(Hy[i]);
    Hr[2][i] = real(Hz[i]);
    Hi[2][i] = imag(Hz[i]);
    
  }

  freeCastMatlab2DArray(Er);
  freeCastMatlab2DArray(Ei);
  freeCastMatlab2DArray(Hr);
  freeCastMatlab2DArray(Hi);
  freeCastMatlab2DArray(vertices);
  free(Ex);
  free(Ey);
  free(Ez);
  free(Hx);
  free(Hy);
  free(Hz);
  free(rad_vec);
  free(theta_vec);
  free(phi_vec);
  free(dims);
}
