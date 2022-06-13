%
%	function [kern,kvals] = plotgridkbkernel(convwidth, overgridfactor, nk)
%
%	Function calculates the actual kernel that is used in the
%	gridding function by repeated calls to gridkb, and plots it.  
%
%	INPUT:
%		convwidth - convwidth to be passed to gridlut.
% 		overgridfactor  Factor by which grid parameters overgrid.
%				(This is gridsize/(FOV/resolution) )
%		nk - number of points at which to evaluate kernel.
%
%	OUTPUT:
%		kern - 2D array of kernel values.
%		kvals - k-space locations corresponding to kern (grid units)
%

%	B. Hargreaves.


function [kern,kvals] = plotlutkernel(convwidth,overgridfactor,nk)


if (nargin < 3)
	nk = 40;
end;
nk = round(nk/2)*2+1;

kvals = 2*(([0:nk-1]/(nk-1))-.5)*convwidth;	% k-points to calculate kernel.
						% units are units of k-space.
kern = zeros(nk,nk);

datsize = 2*ceil(2*convwidth+2+overgridfactor);	% More grid points than width.


% ---- Go through every kernel location, and put a sample
% ---- at that location.  Call gridding function, and keep the
% ---- point at 0,0 as that kernel value.

for m=1:nk
  for n=1:nk
    kloc = (kvals(m)+i*kvals(n)) / (datsize-1) + (eps+i*eps);
    dat = gridkb(kloc,1+i,1,datsize,convwidth,overgridfactor);
    kern(m,n)=dat(datsize/2+1,datsize/2+1);
  end;
end;
    

kern = kern / sqrt(2);		% because 1+i was used as data.
mx = max(abs(kern(:)));

subplot(1,2,1);
c = [1:256]' * [1 1 1]/256;
colormap(c);
image(kvals,kvals,255*abs(kern)/mx);
axis('square');
title('Convolution Kernel from Gridding');
xlabel('Kernel x-location (grid point units)');
ylabel('Kernel y-location (grid point units)');


subplot(1,2,2);
k1d = real(kern((nk+1)/2,:));
plot(real(kvals),k1d);
title('Kernel, along x');
xlabel('k-space location - grid points');


