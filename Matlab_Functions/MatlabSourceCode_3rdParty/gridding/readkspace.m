%
%	function [k,dcf] = readkspace(fname)
%
%	function reads a k-space file of a common format that
%	that is a sequence of float32's with
%	kx1,ky1,dcf1,kx2,ky2,dcf2, ... 
%
%	INPUT:
%		fname	- file name for k-space file.
%
%	OUTPUT
%		k	- k-space points, kx+i*ky.
%		dcf	- density compensation factors.
%

% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: readkspace.m,v $
%	Revision 1.3  2003/09/04 22:10:28  brian
%	Small changes.
%	
%	Revision 1.2  2003/09/04 18:26:01  brian
%	Added warning for k-space < 0.5.
%	
%	Revision 1.1  2002/04/27 20:10:10  bah
%	Added to CVS.
%	
%
% ===========================================================



%

function [k,dcf] = readkspace(fname)

if (nargin < 1)
	fname = 'kspace';
end;

fid = fopen(fname,'r','b');

[arr,n] = fread(fid,inf,'float32');

l = floor(n/3);
arr = arr(1:3*l);

if (l < 1)
	disp('No data in file');
else
	tt = sprintf('%d k-space samples in k-space file %s.',l,fname);
	disp(tt);
end;

arr = reshape(arr,3,l);

k = arr(1,:)+i*arr(2,:);
dcf = arr(3,:);


kmax = max(abs(k));
if (kmax > 0.5)
	disp('Warning:  Maximum k-space > 0.5.  ');
	disp('		Conventionally k-space files should not exceed 0.5');
	disp(' 		This may result in problems later.');
	disp(' ');
	q = input('Normalize k-space to 0.5 max magnitude? (y/n)','s');
	if (q(1)=='y')
		disp('Normalizing...');
		k = 0.48*k/kmax;
	end;
end;

fclose(fid);





