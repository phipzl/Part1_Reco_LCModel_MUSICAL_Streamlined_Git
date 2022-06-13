/*	Gridding Routines.  To be included as header files. */
// 
// 
// #include <math.h>
// #include <stdio.h>
// 
// 
// int imax(x1,x2)
// 
// /* 	Returns maximum of two numbers.
// */
// int x1;
// int x2;
// {
// int retval;
// 
// retval = x1;
// if (x2>x1)
// 	retval = x2;
// 
// printf("Max is %d\n",retval);
// return retval;
// }


// 
// double convkernalpyr(u,v,gridspace)
// 
// /*      Evaluates Convolution kernal at (u,v).  
// 
// 	Pyramid Kernal.
// */
// 
// double u;
// double v;
// double gridspace;
// 
// {
// double retval = 0;
// 
// u = fabs(u)/gridspace;
// v = fabs(v)/gridspace;
// 
// 
// if ((u < 1.0) && (v< 1.0))
//         retval = (1-u) * (1-v);
// else
//         retval = 0.0;
// 
// return retval;
// }


// 
// double convkernalgauss(u,v,width)
// 	
// /*	Evaluates Convolution kernal at (u,v).	*/
// 
// double u;
// double v;
// double width;
// 
// {
// double retval = 0;
// double d;
// double sig;
// 
// d = (u*u+v*v)/width/width;
// 
// if (d < 1.0)
// 	{
// 	sig = (0.0969 + 0.0707*width);
// 	retval = exp(-d/sig);
// 	}
// else
// 	retval = 0.0;
// 
// return retval; 
// }
// 


// 
//  
// double convkernal(u,v,width)
//          
// /*      Evaluates Convolution kernal at (u,v).  */
//  
// double u;
// double v;
// double width;
//  
// {
// double retval = 0;
// double d;
// double sig;
//  
// d = (u*u+v*v)/width/width;
// 
//  
// if (d < 1.0)
//         retval = (1-sqrt(d));
// else
//         retval = 0.0;
// 
// return retval; 
// }
//  





// 
// double ind2k(ind,gridsize)
// 
// /* 	Function converts index value to k-space location. */
// 
// int gridsize;
// int ind;
// 
// {
// double k;
// 
// k = (double)(ind-gridsize/2) / (double)(gridsize);
// return k;
// }
// 
// 
// double k2ind(k,gridsize)
// 
// /*	Function calculates the index of the grid point corresponding
// 	to (k).  Index is fractional - integer index points can
// 	be calculated from this.	*/
// 
// int gridsize;
// double k;
// 
// {
// double gridspace;
// double index;
// 
// gridspace = (double)1.0 / (double)(gridsize);
// index = k/gridspace + (double)(gridsize/2); 
// 
// return index;
// }



// 
// int k2i(k,gridsize)
//  
// /*      Function calculates the index of the grid point corresponding
//         to (k).  Index is fractional - integer index points can
//         be calculated from this.        
// 
// 	Limits to 0:gridsize-1
// */
//  
// int gridsize;
// double k;
//  
// {
// double gridspace;
// int index;
//  
// gridspace = (double)1.0 / (double)(gridsize);
// index = (int)(k/gridspace + (double)(gridsize/2));
// if (index >= gridsize)	
// 	index = gridsize-1;
// else
// 	if (index < 0)
// 		index = 0;
//  
// return index;
// }
//  


// 
// 
// void getpoint(xlow,xhigh,ylow,yhigh,gridsize,sg_real,sg_imag, 
// 		kxx,kyy,rpart,ipart, gridspace)
// 
// /*	Function gets the grid point by interpolating. 	*/
// 
// int xlow;
// int xhigh;
// int ylow;
// int yhigh;
// int gridsize;
// double *sg_real;
// double *sg_imag;
// double kxx;
// double kyy;
// double *rpart;
// double *ipart;
// double gridspace;
// 
// {
// 
// 
// *rpart = 0.0;
// *ipart = 0.0;
// 
// if ((xlow < gridsize) && (xhigh >=0) && (ylow < gridsize) && (yhigh >= 0))
// 	{
// 	if (xlow >= 0)
// 		{
// 		if (ylow >= 0)
// 			{
// 		 	(*rpart) += sg_real[xlow*gridsize+ylow] * convkernal(kxx-ind2k(xlow,gridsize), kyy-ind2k(ylow,gridsize), gridspace); 
// 		 	(*ipart) += sg_imag[xlow*gridsize+ylow] * convkernal(kxx-ind2k(xlow,gridsize), kyy-ind2k(ylow,gridsize), gridspace); 
// 			}
// 		if (yhigh < gridsize)
// 			{
// 		 	(*rpart) += sg_real[xlow*gridsize+yhigh] * convkernal(kxx-ind2k(xlow,gridsize), kyy-ind2k(yhigh,gridsize), gridspace); 
// 		 	(*ipart) += sg_imag[xlow*gridsize+yhigh] * convkernal(kxx-ind2k(xlow,gridsize), kyy-ind2k(yhigh,gridsize), gridspace); 
// 			}
// 		}
// 
// 	if (xhigh < gridsize)
// 		{
// 		if (ylow >= 0)
// 			{
// 		 	(*rpart) += sg_real[xhigh*gridsize+ylow] * convkernal(kxx-ind2k(xhigh,gridsize), kyy-ind2k(ylow,gridsize), gridspace); 
// 		 	(*ipart) += sg_imag[xhigh*gridsize+ylow] * convkernal(kxx-ind2k(xhigh,gridsize), kyy-ind2k(ylow,gridsize), gridspace); 
// 			}
// 		if (yhigh < gridsize)
// 			{
// 		 	(*rpart) += sg_real[xhigh*gridsize+yhigh] * convkernal(kxx-ind2k(xhigh,gridsize), kyy-ind2k(yhigh,gridsize), gridspace); 
// 		 	(*ipart) += sg_imag[xhigh*gridsize+yhigh] * convkernal(kxx-ind2k(xhigh,gridsize), kyy-ind2k(yhigh,gridsize), gridspace); 
// 			}
// 		}
// 
// 	}
// }
// 


// 
// void invgrid(sg_real,sg_imag,gridsize, kx,ky,s_real,s_imag,nsamples)
// 
// 
// double *sg_real;
// double *sg_imag;	 
// int gridsize;		/* Number of points in kx and ky in grid. */
// double *kx;
// double *ky;
// double *s_real;
// double *s_imag;
// int nsamples;		/* Number of samples, total. */
// 
// 
// 
// {
// int count1, count2;
// int gcount1, gcount2;
// 
// double gridspace;
// 
// double kxind,kyind;
// 
// int lowkx;
// int lowky;
// int highkx;
// int highky;
// 
// gridspace = (double)(1.0) / (double)(gridsize);
// 
// 
// for (count1 = 0; count1 < nsamples; count1++)
// 	{
// 
// 	kxind = k2ind(kx[count1],gridsize);
// 	kyind = k2ind(ky[count1],gridsize);
// 
// 	lowkx = (int)kxind;
// 	highkx = (int)(kxind+1.0);
// 	lowky = (int)kyind;
// 	highky = (int)(kyind+1.0);
// 
// 
// 	getpoint(lowkx,highkx,lowky,highky,gridsize,sg_real,sg_imag,
// 		kx[count1],ky[count1],&(s_real[count1]),&(s_imag[count1]),
// 		gridspace);
// 
// 	}
// 
// }
// 


// 
// 
// void calcdensity(kx,ky,nsamples,dens,gridsize,convwidth)
// 
// /*	Function calculates the sample density at each point
// 	using (hopefully) a fast algorithm:
// 
// 	1)  First the density at each sample point is set to 1.
// 	
// 	2)  For each sample in the set S1, we count through the
// 	    samples S2 (which are after S1 in the arrays kx and ky).
// 
// 		If S1 and S2 are within the convolution kernal 
// 		width of each other, the convolution kernal is 
// 		calculated for the vector S1-S2.  This result is
// 		added to the densities at both S1 and S2.
// 
// 		I hope that this minimizes the number of times
// 		that the convolution kernal is calculated!
// 
// 	Brian Hargreaves - Jan 99.
// */
// 
// double *kx;	/* Kx location of each sample point. */
// double *ky;	/* Ky location of each sample point. */
// int nsamples;	/* number of sample points. */
// double *dens;		/* OUTPUT - Density at each point. */
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;	/* Convolution Kernal width in grid points. */
// 
// {
// int count1, count2;
// double gridspace;
//  
// double kxmin,kxmax,kymin,kymax;
// double kxx,kyy;
// double kx1,ky1;
// double kwidth;
// 
// double *pdens1, *pdens2, *pkx, *pky, *pkx1, *pky1;
// double denschange;
// 
// gridspace = (double)(1.0) / (double)(gridsize);
// kwidth = gridspace*convwidth;
// 
// 
// /* ========= Calculate Density Function ======== */
// 
// pdens1 = dens; 
// 
// for (count1 = 0; count1 < nsamples; count1++)
// 	*(pdens1++) = 1.0;
// 
// pdens1 = dens; 
// pkx1 = kx;
// pky1 = ky;
// 
// /* 	Cycle through each sample point (count1), calculating the
// 	density contribution from every other point at this point. */
// 
// for (count1 = 0; count1 < nsamples; count1++)
//         {
// 
// 	/* Find limits so that we don't calculate convolution at
// 	   more points than we need to. */
// 
//         kx1 = *(pkx1++);
//         ky1 = *(pky1++);
//         kxmin = kx1-kwidth;
//         kxmax = kx1+kwidth;
//         kymin = ky1-kwidth;
//         kymax = ky1+kwidth;
// 
// 	pkx = pkx1;
// 	pky = pky1;
// 
// 	/* Get pointer to "remaining points" so that their densities
// 	   can be modified at the same time, since the effect of
// 	   density at S1 by S2 is the same as the effect at S2 by S1. */
// 
// 	pdens2 = pdens1+1;
// 
//         for (count2 = count1+1; count2 < nsamples; count2++)
//                 {
//                 kxx = *(pkx++);
//                 kyy = *(pky++);
//                 if ((kxx>kxmin) && (kxx<kxmax) && (kyy>kymin) && (kyy<kymax))
// 			{
// 			denschange = convkernal(kx1-kxx,
//                                         ky1-kyy, kwidth);
// 			*pdens1 += denschange;
// 			*pdens2 += denschange;
// 			}
// 		pdens2++;
//                 }
// 	pdens1++;
//         }
// 
// 
// 
// }
// 
// 
// void calcraddensity(kx,ky,nsperrad, nsamples,dens,gridsize,convwidth)
// 
// /*	Function calculates the sample density at each point
// 	using (hopefully) a fast algorithm:
// 
// 	1)  First the density at each sample point is set to 1.
// 	
// 	2)  For each sample in the set S1, we count through the
// 	    samples S2 (which are after S1 in the arrays kx and ky).
// 
// 		If S1 and S2 are within the convolution kernal 
// 		width of each other, the convolution kernal is 
// 		calculated for the vector S1-S2.  This result is
// 		added to the densities at both S1 and S2.
// 
// 		I hope that this minimizes the number of times
// 		that the convolution kernal is calculated!
// 
// 	Brian Hargreaves - Jan 99.
// */
// 
// double *kx;	/* Kx location of each sample point. */
// double *ky;	/* Ky location of each sample point. */
// int nsperrad;	/* number of sample points per radial trajectory. */
// int nsamples;	/* number of sample points. */
// double *dens;		/* OUTPUT - Density at each point. */
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;	/* Convolution Kernal width in grid points. */
// 
// {
// int count1, count2;
// double gridspace;
//  
// double kxmin,kxmax,kymin,kymax;
// double kxx,kyy;
// double kx1,ky1;
// double kwidth;
// 
// double *pdens1, *pdens2, *pkx, *pky, *pkx1, *pky1;
// double denschange;
// int ntheta;
// 
// gridspace = (double)(1.0) / (double)(gridsize);
// kwidth = gridspace*convwidth;
// 
// 
// /* ========= Calculate Density Function ======== */
// 
// pdens1 = dens; 
// 
// for (count1 = 0; count1 < nsamples; count1++)
// 	*(pdens1++) = 1.0;
// 
// pdens1 = dens; 
// pkx1 = kx;
// pky1 = ky;
// 
// /* 	Cycle through each sample point (count1), calculating the
// 	density contribution from every other point at this point. */
// 
// for (count1 = 0; count1 < nsperrad; count1++)
//         {
// 
// 	/* Find limits so that we don't calculate convolution at
// 	   more points than we need to. */
// 
//         kx1 = *(pkx1++);
//         ky1 = *(pky1++);
//         kxmin = kx1-kwidth;
//         kxmax = kx1+kwidth;
//         kymin = ky1-kwidth;
//         kymax = ky1+kwidth;
// 
// 	pkx = pkx1;
// 	pky = pky1;
// 
// 	/* Get pointer to "remaining points" so that their densities
// 	   can be modified at the same time, since the effect of
// 	   density at S1 by S2 is the same as the effect at S2 by S1. */
// 
// 	pdens2 = pdens1+1;
// 
//         for (count2 = count1+1; count2 < nsamples; count2++)
//                 {
//                 kxx = *(pkx++);
//                 kyy = *(pky++);
//                 if ((kxx>kxmin) && (kxx<kxmax) && (kyy>kymin) && (kyy<kymax))
// 			{
// 			denschange = convkernal(kx1-kxx,
//                                         ky1-kyy, kwidth);
// 			*pdens1 += denschange;
// 			*pdens2 += denschange;
// 			}
// 		pdens2++;
//                 }
// 	pdens1++;
//         }
// 
// /*	Now duplicate the densities for other angles... */
// 
// ntheta = nsamples/nsperrad;
// 
// 	
// 	
// pdens2 = dens+nsperrad; 
// 	
// for (count1 = 1; count1 < ntheta; count1++)
// 	{
// 	pdens1 = dens; 
// 	for (count2 = 0; count2 < nsperrad; count2++)
// 		*(pdens2++) = *(pdens1++);
// 	}
// 
// }
// 

// 
// void gridonly(kx,ky,s_real,s_imag,nsamples, dens,
// 	sg_real,sg_imag, gridsize, convwidth)
// 
// double *kx;
// double *ky;
// double *s_real;		/* Sampled data. */
// double *s_imag;
// int nsamples;		/* Number of samples, total. */
// 
// double *dens;		/* Output - density function. */
// 
// double *sg_real;
// double *sg_imag;	 
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;
// 
// 
// {
// int count1, count2;
// int gcount1, gcount2;
// double gx,gy;
// int col;
// int loc;
// double gridspace;
// 
// double kxmin,kxmax,kymin,kymax;
// double kxx,kyy;
// double kx1,ky1;
// double kwidth;
// int ixmin,ixmax,iymin,iymax;
// 
// 
// gridspace = (double)(1.0) / (double)(gridsize);
// kwidth = gridspace*convwidth;
// 
// /* ========= Zero Output Points ========== */
//  
// for (gcount1 = 0; gcount1 < gridsize; gcount1++)
//         {
//         col = gcount1*gridsize;
//         for (gcount2 = 0; gcount2 < gridsize; gcount2++)
//                 {
//                 loc = col+gcount2;
//                 sg_real[loc] = 0.0;
//       		sg_imag[loc] = 0.0;
// 		}
// 	} 
//  
// /* ========= Find Grid Points ========= */
// 
// 
//                 
// for (count1 = 0; count1 < nsamples; count1++)
// 	{
// 		/* Box around k sample location. */
// 
// 	kxx = kx[count1];
// 	kyy = ky[count1];
// 	ixmin = k2i(kxx - kwidth,gridsize);
// 	ixmax = k2i(kxx + kwidth,gridsize);
// 	iymin = k2i(kyy - kwidth,gridsize);
// 	iymax = k2i(kyy + kwidth,gridsize);
// 
// 	/*
// 	printf("min(%d,%d) - max(%d,%d) \n",ixmin,iymin,ixmax,iymax);
// 	*/
// 
// 	for (gcount1 = ixmin; gcount1 <= ixmax; gcount1++)
//         	{
// 		gx = (double)(gcount1-gridsize/2) / (double)gridsize;
//         	col = gcount1*gridsize;
// 		for (gcount2 = iymin; gcount2 <= iymax; gcount2++)
// 			{
//                 	gy = (double)(gcount2-gridsize/2) / (double)gridsize;
// 			/*
// 			printf("Separation (%f,%f), c(k)= %f\n",gx-kxx,gy-kyy,
// 				convkernal(gx-kxx, gy-kyy,kwidth));
// 			*/
// 
// 			sg_real[col+gcount2] += convkernal(gx-kxx, gy-kyy,
//                                 kwidth) * s_real[count1]/dens[count1];
// 			sg_imag[col+gcount2] += convkernal(gx-kxx, gy-kyy,
//                                 kwidth) * s_imag[count1]/dens[count1];
// 			}
// 		}
// 	}
// }
// 
// 
// 
// void calcdcflut(kx,ky,nsamples,dcf,gridsize,convwidth,kerneltable,nkernelpts)
// 
// /* Calculation of density compensation factors using 
// 			lookup-table convolution kernel */
// 	/* See Notes below */
// 
// /* INPUT/OUTPUT */
// 
// double *kx;		/* Array of kx locations of samples. */
// double *ky;		/* Array of ky locations of samples. */
// int nsamples;		/* Number of k-space samples, total. */
// double *dcf;		/* Output: Density compensation factors. */
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;	/* Kernel width, in grid points.	*/
// double *kerneltable;	/* 1D array of convolution kernel values, starting
// 				at 0, and going to convwidth. */
// int nkernelpts;		/* Number of points in kernel lookup-table */
// 
// /*------------------------------------------------------------------
// 	NOTES:
// 
// 	Assume gridding uses the following formula, which 
// 	describes the contribution of each data sample 
// 	to the value at each grid point:
// 
// 		grid-point value += data value / dcf * kernel(dk)
// 
// 	where:
// 		data value is the complex sampled data value.
// 		dcf is the density compensation factor for the sample point.
// 		kernel is the convolution kernel function.
// 		dk is the k-space separation between sample point and
// 			grid point.
// 
// 	"grid units"  are integers from 0 to gridsize-1, where
// 	the grid represents a k-space of -.5 to .5.
// 
// 	The convolution kernel is passed as a series of values 
// 	that correspond to "kernel indices" 0 to nkernelpoints-1.
//   ------------------------------------------------------------------ */
// 
// {
// int kcount1;		/* Counts through k-space sample locations */
// int kcount2;		/* Counts through k-space sample locations */
// 
// double *dcfptr;			/* Aux. pointer, for looping.	*/
// double *kxptr;			/* Aux. pointer. */
// double *kyptr;			/* Aux. pointer. */
// 
// double kwidth;			/* Conv kernel width, in k-space units */
// double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
// int kernelind;			/* Index of kernel value, for lookup. */
// double fracdk;			/* Fractional part of lookup index. */
// double fracind;			/* Fractional lookup index. */
// double kern;			/* Kernel value, avoid duplicate calculation.*/
// 
// kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */
// 
// 
// /* ========= Set all DCFs to 1. ========== */
// 
// /* 	DCF = 1/(sum( ksamps convolved with kernel ) 	*/
// 
// dcfptr = dcf; 
// for (kcount1 = 0; kcount1 < nsamples; kcount1++)
//         *dcfptr++ = 1.0;
// 
//  
// /* ========= Loop Through k-space Samples ========= */
//                 
// dcfptr = dcf;
// kxptr = kx;
// kyptr = ky;
//  
// for (kcount1 = 0; kcount1 < nsamples; kcount1++)
// 	{
// 	/* printf("Current k-space location = (%f,%f) \n",*kxptr,*kyptr); */
// 
// 	for (kcount2 = kcount1+1; kcount2 < nsamples; kcount2++)
// 		{
// 		dkx = *kxptr-kx[kcount2];
// 		dky = *kyptr-ky[kcount2];
// 		dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/
// 		/*
// 		printf("   Comparing with k=(%f,%f),  sep=%f \n",
// 				kx[kcount2],ky[kcount2],dk);	
// 		*/
// 	
// 		if (dk < kwidth)	/* sample affects this
// 					   grid point		*/
// 		    {
// 			/* Find index in kernel lookup table */
// 		    fracind = dk/kwidth*(double)(nkernelpts-1);
// 		    kernelind = (int)fracind;
// 		    fracdk = fracind-(double)kernelind;
// 
// 			/* Linearly interpolate in kernel lut */
// 		    kern = kerneltable[(int)kernelind]*(1-fracdk)+
// 		    		kerneltable[(int)kernelind+1]*fracdk;
// 
// 		    dcf[kcount2] += kern;
// 		    *dcfptr += kern;
// 		    }
// 		}
// 	dcfptr++;
// 	kxptr++;
// 	kyptr++;
// 	}
// 
// 
// dcfptr = dcf; 
// for (kcount1 = 0; kcount1 < nsamples; kcount1++)
//         *dcfptr = 1.0/ *dcfptr++;
// 
// }


// 
// void gridlut(kx,ky,s_real,s_imag,nsamples, dcf,
// 	sg_real,sg_imag, gridsize, convwidth, kerneltable, nkernelpts)
// 
// /* Gridding function that uses a lookup table for a circularyl
// 	symmetric convolution kernel, with linear interpolation.  
// 	See Notes below */
// 
// /* INPUT/OUTPUT */
// 
// double *kx;		/* Array of kx locations of samples. */
// double *ky;		/* Array of ky locations of samples. */
// double *s_real;		/* Sampled data, real part. */
// double *s_imag; 	/* Sampled data, real part. */
// int nsamples;		/* Number of k-space samples, total. */
// double *dcf;		/* Density compensation factors. */
// 
// double *sg_real;	/* OUTPUT array, real parts of data. */
// double *sg_imag;	/* OUTPUT array, imag parts of data. */	 
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;	/* Kernel width, in grid points.	*/
// double *kerneltable;	/* 1D array of convolution kernel values, starting
// 				at 0, and going to convwidth. */
// int nkernelpts;		/* Number of points in kernel lookup-table */
// 
// /*------------------------------------------------------------------
// 	NOTES:
// 
// 	This uses the following formula, which describes the contribution
// 	of each data sample to the value at each grid point:
// 
// 		grid-point value += data value / dcf * kernel(dk)
// 
// 	where:
// 		data value is the complex sampled data value.
// 		dcf is the density compensation factor for the sample point.
// 		kernel is the convolution kernel function.
// 		dk is the k-space separation between sample point and
// 			grid point.
// 
// 	"grid units"  are integers from 0 to gridsize-1, where
// 	the grid represents a k-space of -.5 to .5.
// 
// 	The convolution kernel is passed as a series of values 
// 	that correspond to "kernel indices" 0 to nkernelpoints-1.
//   ------------------------------------------------------------------ */
// 
// {
// int kcount;		/* Counts through k-space sample locations */
// int gcount1, gcount2;	/* Counters for loops */
// int col;		/* Grid Columns, for faster lookup. */
// 
// double kwidth;			/* Conv kernel width, in k-space units */
// double dkx,dky,dk;		/* Delta in x, y and abs(k) for kernel calc.*/
// int ixmin,ixmax,iymin,iymax;	/* min/max indices that current k may affect*/
// int kernelind;			/* Index of kernel value, for lookup. */
// double fracdk;			/* Fractional part of lookup index. */
// double fracind;			/* Fractional lookup index. */
// double kern;			/* Kernel value, avoid duplicate calculation.*/
// double *sgrptr;			/* Aux. pointer, for loop. */
// double *sgiptr;			/* Aux. pointer, for loop. */
// int gridsizesq;			/* Square of gridsize */
// int gridcenter;			/* Index in output of kx,ky=0 points. */
// int gptr_cinc;			/* Increment for grid pointer. */
// gridcenter = gridsize/2;	/* Index of center of grid. */
// kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */
// 
// 
// /* ========= Zero Output Points ========== */
// 
// sgrptr = sg_real;
// sgiptr = sg_imag;
// gridsizesq = gridsize*gridsize;
//  
// for (gcount1 = 0; gcount1 < gridsizesq; gcount1++)
//         {
// 	*sgrptr++ = 0;
// 	*sgiptr++ = 0;
// 	}
// 
//  
// /* ========= Loop Through k-space Samples ========= */
//                 
// for (kcount = 0; kcount < nsamples; kcount++)
// 	{
// 
// 	/* ----- Find limit indices of grid points that current
// 		 sample may affect (and check they are within grid) ----- */
// 
// 	ixmin = (int) ((*kx-kwidth)*gridsize +gridcenter);
// 	if (ixmin < 0) ixmin=0;
// 	ixmax = (int) ((*kx+kwidth)*gridsize +gridcenter)+1;
// 	if (ixmax >= gridsize) ixmax=gridsize-1;
// 	iymin = (int) ((*ky-kwidth)*gridsize +gridcenter);
// 	if (iymin < 0) iymin=0;
// 	iymax = (int) ((*ky+kwidth)*gridsize +gridcenter)+1;
// 	if (iymax >= gridsize) iymax=gridsize-1;
// 
// 		
// 	  /* Increment for grid pointer at end of column to top of next col.*/
// 	gptr_cinc = gridsize-(iymax-iymin)-1;	/* 1 b/c at least 1 sgrptr++ */
// 		
// 	sgrptr = sg_real + (ixmin*gridsize+iymin);
// 	sgiptr = sg_imag + (ixmin*gridsize+iymin);
// 						
// 	for (gcount1 = ixmin; gcount1 <= ixmax; gcount1++)
//         	{
// 		dkx = (double)(gcount1-gridcenter) / (double)gridsize  - *kx;
// 		for (gcount2 = iymin; gcount2 <= iymax; gcount2++)
// 			{
//                 	dky = (double)(gcount2-gridcenter) / 
// 					(double)gridsize - *ky;
// 
// 			dk = sqrt(dkx*dkx+dky*dky);	/* k-space separation*/
// 
// 			if (dk < kwidth)	/* sample affects this
// 						   grid point		*/
// 			    {
// 				/* Find index in kernel lookup table */
// 			    fracind = dk/kwidth*(double)(nkernelpts-1);
// 			    kernelind = (int)fracind;
// 			    fracdk = fracind-(double)kernelind;
// 
// 				/* Linearly interpolate in kernel lut */
// 			    kern = kerneltable[(int)kernelind]*(1-fracdk)+
// 			    		kerneltable[(int)kernelind+1]*fracdk;
// 
// 			    *sgrptr += kern * *s_real * *dcf;
// 			    *sgiptr += kern * *s_imag * *dcf;
// 			    }
// 			sgrptr++;
// 			sgiptr++;
// 			}
// 		sgrptr+= gptr_cinc;
// 		sgiptr+= gptr_cinc;
// 		}
// 	kx++;		/* Advance kx pointer */
// 	ky++;   	/* Advance ky pointer */
// 	dcf++;		/* Advance dcf pointer */		
// 	s_real++;	/* Advance real-sample pointer */
// 	s_imag++;	/* Advance imag-sample pointer */
// 	}
// }
// 


// 
// void griddata(kx,ky,s_real,s_imag,nsamples, dens,
// 	sg_real,sg_imag, gridsize, convwidth)
// 
// double *kx;
// double *ky;
// double *s_real;
// double *s_imag;
// int nsamples;		/* Number of samples, total. */
// 
// double *dens;		/* Output - density function. */
// 
// double *sg_real;
// double *sg_imag;	 
// int gridsize;		/* Number of points in kx and ky in grid. */
// double convwidth;
// 
// 
// {
// int count1, count2;
// int gcount1, gcount2;
// double gx,gy;
// int col;
// int loc;
// double gridspace;
// 
// double kxmin,kxmax,kymin,kymax;
// double kxx,kyy;
// double kx1,ky1;
// double kwidth;
// 
// gridspace = (double)(1.0) / (double)(gridsize);
// kwidth = gridspace*convwidth;
// 
// 
// /* ========= Calculate Density Function ======== */
// 
// calcdensity(kx,ky,nsamples,dens,gridsize,convwidth);
// 
// 
// 
// /* ========= Find Grid Points ========= */
// 
// for (gcount1 = 0; gcount1 < gridsize; gcount1++)
// 	{
// 	gx = (double)(gcount1-gridsize/2) / (double)gridsize;
// 	col = gcount1*gridsize;
// 		
// 	kxmin = gx-kwidth;
// 	kxmax = gx+kwidth;
// 
// 	for (gcount2 = 0; gcount2 < gridsize; gcount2++)
// 		{
// 		loc = col+gcount2;
// 		sg_real[loc] = 0.0;
// 		sg_imag[loc] = 0.0;
// 
// 		gy = (double)(gcount2-gridsize/2) / (double)gridsize;
// 
// 		kymin = gy-kwidth;
// 		kymax = gy+kwidth;
// 
// 		for (count1 = 0; count1 < nsamples; count1++)
// 			{
// 			kxx = kx[count1];
// 			kyy = ky[count1];
// 
// 			if ((kxx>kxmin)&&(kxx<kxmax)&&(kyy>kymin)&&(kyy<kymax))
// 			    {
// 			    sg_real[loc] += convkernal(gx-kxx, gy-kyy,
// 				gridspace*convwidth) * s_real[count1]/dens[count1];
// 			    sg_imag[loc] += convkernal(gx-kxx, gy-kyy,
// 				gridspace*convwidth) * s_imag[count1]/dens[count1];
// 			    }
// 			}
// 		}
// 
// 	}
// 
// }
// 
// 
// 

