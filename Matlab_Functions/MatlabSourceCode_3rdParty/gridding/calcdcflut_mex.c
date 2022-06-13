// 
// #include "mex.h" 
// #include "gridroutines.c"
// 
// 
// 
// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// 
// lllllllllllllllllllllllll
// {
// double *kx;
// double *ky;
// int nsamples;
// double *dcf;
// int gridsize;
// double convwidth;
// double *kerneltable;
// int nkernelpts;
// 
// 
// nsamples = mxGetM(prhs[0]) * mxGetN(prhs[0]);
// 
// kx = mxGetPr(prhs[0]);
// ky = mxGetPi(prhs[0]);
// gridsize = (int)(*mxGetPr(prhs[1]));
// 	
// if (nrhs > 3)
// 	convwidth = *mxGetPr(prhs[2]);
// else
// 	convwidth = 1.0;
// 
// kerneltable = mxGetPr(prhs[3]);
// nkernelpts = mxGetM(prhs[3]) * mxGetN(prhs[3]);
// 
// plhs[0] = mxCreateDoubleMatrix(nsamples,1,mxREAL);
// dcf = mxGetPr(plhs[0]);
// 
// calcdcflut(kx,ky,nsamples,dcf,gridsize,convwidth,kerneltable,nkernelpts);
// 
// }
// 
// 



