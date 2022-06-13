%
%   [dat,kg]  = plotgridkb(ktraj,dcf,gridsize,kwidth,overgridfactor,K)
%
%	Function grids 1's onto the given trajectory using the
%	appropriate Kaiser-Bessel kernel given
%	the convolution width and the overgridding factor.
%
%	INPUT:
%		ktraj	- k-space trajectory.
%		dcf	- density compensation factors.
%		gridsize- size of grid to (normally) grid to.
%		kwidth  - convolution kernel width.
%		overgridfactor - alpha parameter for gridding.
%		K 	- #points between grid points for plot.
%
%	OUTPUT:
%		kg	- corresponding k-space for calculated locations.
%		dat	- gridded data.
%
%
%	B.Hargreaves.



function [dat,kg]  = plotgridkb(ktraj,dcf,gridsize,kwidth,overgridfactor,K)

ksamps = ones(length(ktraj(:))) * [1+eps*i];

if (nargin < 6)
	K = 8;		% oversampling for plotting.
end;

% ========== Calculate Kaiser-Bessel Kernel for given parameters =========

kernellength = 32;					% Arbitrary.
u = [0:kernellength-1]/(kernellength-1) * kwidth/2;	% Kernel radii.
kerneltable = 0*u;					% Allocate.

kerneltable = kb(u,kwidth, pi*kwidth*(overgridfactor-.5) );
%for k=1:length(u)
%	kerneltable(k) = kb(u(k),kwidth, pi*kwidth*(overgridfactor-.5) );
%end;

kerneltable=calckbkernel(kwidth,overgridfactor);	% Calculate kernel.



% ========== Do Gridding =========

dat = gridlut(ktraj,ksamps,dcf,K*gridsize,K*kwidth/2, kerneltable);

dat = dat/max(abs(dat(:)));


% ========== Highlight Grid Locations =====

for m=0:gridsize-1
  for n=0:gridsize-1
    dat(m*K+1,n*K+1) = 1-round(abs(dat(m*K+1,n*K+1)));
  end;
end;

kg = [0:K*gridsize-1]/gridsize/K-.5;

% ========= Plot grid. =========

c = [1:256]' * [1 1 1]/256;
colormap(c);
image(kg,kg,256*abs(dat));
title('Gidding of Ones at Sample Points');
xlabel('k_x (in units of k-space)');
ylabel('k_y (in units of k-space)');
 



