/*************************************************************************************
	
				  Implements the Richards-Wolf integrals

					Development : P. Munro 

**************************************************************************************/
#include "math.h"
#include "mex.h"
#include <complex>
#include "matrix.h"
using namespace std;
#include "globals.h"
#include "matlabio.h"



#ifdef __linux
#else
#define M_PI 3.14159265358979
#endif

// ********************************************************
//			Impelmentation of functions
// ********************************************************


void gauleg(double x1, double x2, double x[], double w[], int n);

/*Calculate the I values for RW diffraction integral*/

struct complex_vector Richards_Wolf(double r, double z, double a1, double a2, double k, double F, int n) 
{ 

        #define                         NIMAX 10000 
		//defined somewhere as 300, could use dynamic but this is quicker.
        double                  abscissa[NIMAX]; 
        double                  weight[NIMAX]; 

        struct complex_vector ret; 
        gauleg(a1,a2,abscissa,weight,n); 
					
        complex<double> I0 = 0.0; 
        complex<double> I1 = 0.0; 
        complex<double> I2 = 0.0; 

        for (int i=0;i<n;i++) { 

                // definition of angles 

                double st = sin(abscissa[i]); 
                double ct = cos(abscissa[i]); 

                // amplitude terms 

                double ai0 = sqrt(ct)*(1.0+ct)*st; 
                double ai1 = sqrt(ct)*st*st; 
                double ai2 = sqrt(ct)*(1.0-ct)*st; 

                // arguments of functions 

                double argj = k*r*st; 
                double arge = k*z*ct; 

                // functions 

                complex<double> ctmp = 2.*M_PI*exp(I*arge)*exp(-pow(st*F/sin(a2),2))*weight[i]; 

                I0 +=ai0*j0(argj)*ctmp; 
                I1 +=ai1*j1(argj)*ctmp; 
                I2 +=ai2*jn(2,argj)*ctmp; 
        } 

        ret.X = I0; 
        ret.Y = I1; 
        ret.Z = I2; 

        return(ret); 

} 

#define EPS 3.0e-13

/*Set up the integration routine parameters*/

void gauleg(double x1, double x2, double x[], double w[], int n)
{
	int		m, j, i;
	double	z1 ,z , xm, xl, pp, p3, p2, p1;

	m  = (n+1) / 2;
	xm = 0.5 * (x2+x1);
	xl = 0.5 * (x2-x1);

	for (i = 1; i <= m; i++)  
	{
		z = cos( dcpi * (i-0.25) / (n+0.5) );
		do 
		{
			p1 = 1.0;
			p2 = 0.0;
			for ( j = 1; j <= n; j++ )
			{
				p3 = p2;
				p2 = p1;
				p1 = ( (2.0*j-1.0)*z*p2- (j-1.0)*p3 ) / j;
			}
			pp = n * (z*p1-p2) / (z*z-1.0);
			z1 = z;
			z  = z1 - p1/pp;
		} 
		while ( fabs(z-z1) > EPS );
		x[i-1] = xm - xl*z;
		x[n-i] = xm + xl*z;
		w[i-1] = 2.0*xl / ((1.0-z*z)*pp*pp);
		w[n-i] = w[i-1];
	}
}

#undef EPS


/*The interface to matlab. This function is called by matlab and it does the necessary data structure conversion 
  and creation*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
	/*Expect 6 inputs:

	  NA - numerical aperture of focusing lens
	  k  - wavenumber
	  xverts - x component of vertices
	  yverts - y component of vertices
	  zverts - z component of vertices
	  F - parameter of Gaussian apodisation
	  nintegral - number of integration steps
	  
	  Expect 3 outputs:

	  Ex
	  Ey 
	  Ez

	*/
  double NA, knum, *xverts, *yverts, *zverts, *Exr, *Exi, *Eyr, *Eyi, *Ezr, *Ezi, x, y, z, phi, F;
  int nintegral, count = 0, nverts, ndim;
  const mwSize *dims;
  struct complex_vector ret; 
  complex<double> I0, I1, I2, Ex, Ey, Ez;
	
	//first check that we have 11 inputs
	if (nrhs != 7) 
		mexErrMsgTxt("7 inputs required.");
	//and number of outputs
	if (nlhs != 3) 
	    mexErrMsgTxt("3 outputs required.");

	NA = *mxGetPr(prhs[count++]);
	knum = *mxGetPr(prhs[count++]);
	xverts = mxGetPr(prhs[count]); 

	ndim = mxGetNumberOfDimensions(prhs[count]); 
	dims = mxGetDimensions(prhs[count]); 

	nverts = mxGetNumberOfElements(prhs[count++]);
	yverts = mxGetPr(prhs[count]);
	if(nverts != mxGetNumberOfElements(prhs[count++]))
	  mexErrMsgTxt("y vertice list has incorrect dimension");
	zverts = mxGetPr(prhs[count]);
	if(nverts != mxGetNumberOfElements(prhs[count++]))
	  mexErrMsgTxt("z vertice list has incorrect dimension");
	F = (double)*mxGetPr(prhs[count++]);
	nintegral = (int)*mxGetPr(prhs[count++]);

	fprintf(stdout, "F: %e\n",F);
	if(ReportVerbosely()){
	  fprintf(stdout, "NA: %e\n", NA);
	  fprintf(stdout, "k: %e\n", knum);
	  fprintf(stdout, "nverts: %d\n",nverts);
	  
	  fprintf(stdout, "nintegral: %d\n", nintegral);
	}
	
	plhs[0] = mxCreateNumericArray( ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
	plhs[1] = mxCreateNumericArray( ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
	plhs[2] = mxCreateNumericArray( ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);

	Exr = mxGetPr(plhs[0]);
	Exi = mxGetPi(plhs[0]);
	Eyr = mxGetPr(plhs[1]);
	Eyi = mxGetPi(plhs[1]);
	Ezr = mxGetPr(plhs[2]);
	Ezi = mxGetPi(plhs[2]);

	//now loop over each (x,y) pair
	for(int i = 0; i < nverts; i++){
	  x = xverts[i];
	  y = yverts[i];
	  z = zverts[i];

	  phi = atan2(y,x);

	  ret = Richards_Wolf( sqrt(x*x + y*y), z, 0, asin(NA), knum, F, nintegral);
	  I0 = ret.X;
	  I1 = ret.Y;
	  I2 = ret.Z;

	  Ex = -1.*I*(I0 + I2*cos(2.*phi));
	  Ey = -1.*I*(I2*sin(2.*phi));
	  Ez = -1.*I*(-2.*I*cos(phi)*I1);

	  Exr[i] = real(Ex);
	  Exi[i] = imag(Ex);
	  Eyr[i] = real(Ey);
	  Eyi[i] = imag(Ey);
	  Ezr[i] = real(Ez);
	  Ezi[i] = imag(Ez);
	}
}


/*void main(int n){
	struct complex_vector retval;
	
	retval = Richards_Wolf(.1, 100, 1,1e7 , 100); 
	fprintf(stderr, "\n%f+i%f %f+i%f %f+i%f\n",retval.X.real(),retval.X.imag(),retval.Y.real(),retval.Y.imag(),retval.Z.real(),retval.Z.imag());

}*/

