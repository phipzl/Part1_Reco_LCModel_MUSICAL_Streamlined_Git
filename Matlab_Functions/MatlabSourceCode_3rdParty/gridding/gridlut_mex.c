// 
// #include "mex.h" 
// #include "gridroutines.c"
// 
// 
// 
// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// 
// 
// {
// double *kx;		/* k-space trajectory x locations.	*/
// double *ky;		/* k-space trajectory y locations.	*/
// double *s_real;		/* real-part of data samples.	*/
// double *s_imag;		/* imaginary-part of data samples. */
// int nsamples;		/* number of data samples (or trajectory points) */
// double *dens;		/* (pre) density correction factor at each traj. loc.*/
// double *sg_real;	/* real-part of gridded data sample. */
// double *sg_imag;	/* imaginary-part of gridded data sample. */
// int gridsize;		/* size of grid, in samples.  ie gridsize x gridsize */
// double convwidth;	/* width of convolution kernel for gridding. */
// double *kerneltable;	/* lookup-table for linearly-interpolated kernel.  */
// int nkernelpts;		/* number of points in kernel lookup table. */
// 
// mxArray *dens_mat;	
// 
// nsamples = mxGetM(prhs[0]) * mxGetN(prhs[0]);	/* Samples may be passed as
// 							1xN, Nx1 or 2D array */
// 
// kx = mxGetPr(prhs[0]);		/* Get kx locations */
// ky = mxGetPi(prhs[0]);		/* Get ky locations */
// s_real = mxGetPr(prhs[1]);	/* Get real parts of data samples. */
// s_imag = mxGetPi(prhs[1]);	/* Get imaginary parts of data samples. */
// dens = mxGetPr(prhs[2]);	/* Get density correction factors. */
// gridsize = (int)(*mxGetPr(prhs[3]));	/* Get grid size. */
// 	
// if (nrhs > 4)
// 	convwidth = *mxGetPr(prhs[4]);			/* Get conv width */
// else
// 	convwidth = 1.0;				/* Assign default. */
// 
// kerneltable = mxGetPr(prhs[5]);				/* Get kernel table.*/
// nkernelpts = mxGetM(prhs[5]) * mxGetN(prhs[5]);		/* and # points. */
// 
// plhs[0] = mxCreateDoubleMatrix(gridsize,gridsize,mxCOMPLEX);	/* output */
// sg_real = mxGetPr(plhs[0]);
// sg_imag = mxGetPi(plhs[0]);
// 
// 	
// gridlut(kx,ky,s_real,s_imag,nsamples, dens,
// 		sg_real,sg_imag, gridsize, convwidth, kerneltable, nkernelpts);
// 
// }
// 
// 
// 
// 

