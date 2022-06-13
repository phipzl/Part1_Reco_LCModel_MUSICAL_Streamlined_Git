%
%	function [kern,kvals] = plotlutkernel(convwidth, kerneltable, nk)
%
%	Function calculates the actual kernel that is used in the
%	gridding function by repeated calls to gridlut, and plots it.  
%
%	INPUT:
%		convwidth - convwidth to be passed to gridlut.
%		kerneltable - kerneltable to be passed to gridlut.
%		nk - number of points at which to evaluate kernel.
%
%	OUTPUT:
%		kern - 2D array of kernel values.
%		kvals - k-space locations corresponding to kern (grid units)
%

%	B. Hargreaves.


function [kern,kvals] = plotlutkernel(convwidth,kerneltable,nk)


if (nargin < 3)
	nk = 20;
end;


kvals = 2*(([0:nk-1]/(nk-1))-.5)*convwidth;	% k-points to calculate kernel.
						% units are units of k-space.
kern = zeros(nk,nk);

datsize = ceil(2*convwidth+3);	% More grid points than width.
for m=1:nk
  for n=1:nk
    kloc = (kvals(m)+i*kvals(n)) / (datsize-1) + (eps+i*eps);
    dat = gridlut(kloc,1+i,1,datsize,convwidth,kerneltable);
    kern(m,n)=dat(datsize/2+1,datsize/2+1);
  end;
end;
    

kern = kern / sqrt(2);		% because 1+i was used as data.

mx = max(abs(kern(:)));

c = [1:256]' * [1 1 1]/256;
colormap(c);
image(kvals,kvals,255*abs(kern)/mx);
axis('square');
title('Convolution Kernel from Gridding');
xlabel('Kernel x-location (grid point units)');
ylabel('Kernel y-location (grid point units)');



