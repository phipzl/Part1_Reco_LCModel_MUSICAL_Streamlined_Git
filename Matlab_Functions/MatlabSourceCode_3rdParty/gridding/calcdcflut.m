%
%	function [dcf] = calcdcflut(ktraj,gridsize,convwidth, kerneltable) 
%
%	Function calculates the density compensation factors using
%	the convolution method.
%
%	ktraj		K-space sample locations (complex).
%	gridsize 	Size of grid onto which to grid data.
%	convwidth	(optional) width of convolution kernal in grid points.
%	kerneltable	Convolution kernel points from 0 to convwidth.
%
%

function [dcf] = calcdcflut(ktraj,gridsize,convwidth, kerneltable) 



% --------- Check k-space trajectory is complex-valued -------
s = size(ktraj);
if (s(2)==2)
	ktraj = ktraj(:,1)+i*ktraj(:,2);
end;
if isreal(ktraj)
	ktraj = ktraj + 2*eps*i;
	disp('Warning: k-space trajectory should be complex-valued, or Nx2');
end;



% --------- Check k-space trajectory does not exceed 0.5 -------
if ( max(abs(real(ktraj)))>=0.5)  | (max(abs(imag(ktraj)))>=0.5 )
	disp('Warning:  k-space location radii should not exceed 0.5.');
	disp('		Some samples could be ignored.');
end;



% --------- Default Arguments -------------
if (nargin < 2)
	gridsize=256;
end;
if (nargin < 3)
	convwidth=2;
end;
if (nargin < 4)
	kerneltable=[1 0];	 
end;

% --------- Call Mex Function for this ---------


dcf = calcdcflut_mex(ktraj,gridsize,convwidth,kerneltable);



