%
%	function writekspace(fname,k,dcf)
%
%	function writes a k-space file of the commonly-used 
%	format that is a sequence of float32's with
%	kx1,ky1,dcf1,kx2,ky2,dcf2, ... 
%
%	INPUT:
%		fname	- file name for k-space file.
%		k	- k-space points, kx+i*ky.
%		dcf	- density compensation factors.
%

% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: writekspace.m,v $
%	Revision 1.2  2003/02/11 19:41:55  brian
%	Added a bunch of signal simulation files.
%	Not quite sure what the changes to other
%	files were here.
%	
%	Revision 1.1  2002/04/27 20:10:10  bah
%	Added to CVS.
%	
%
% ===========================================================



function writekspace(fname,k,dcf)

if (max(abs(k(:))) > 0.5)
	disp('Warning:  k-space amplitude exceeds 0.5');
end;

fid = fopen(fname,'w');

ll = length(k);
k=reshape(k,1,ll);
dcf=reshape(dcf,1,ll);

arr = [real(k);imag(k); dcf];

count = fwrite(fid,arr,'float');

if (count < 3*ll)
	disp('Error writing file');
end;

fclose(fid);





