function [x] = GaussSeidelInversion(A, b, x0, iters)%
% ExponentialFilter Apply an Exponential filter to time-domain signals (e.g. FIDs)
%
%
%
%
% The function computes an exponential filter in Hertz
%
%
% [OutArray,exp_filter_funct] = ExponentialFilter(InArray,dwelltime,ApplyAlongDim,exp_filter_Hz)
%
% Input: 
% -         A                     ...    bla
%
% Output:
% -         OutArray                    ...     bla
%
%
% File dependancy:

% Further remarks: 
% This was copied from: http://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method#Program_to_solve_arbitrary_no._of_equations_using_Matlab


%% 0. Declarations, Preparations, Definitions





%% 1.


n = length(A);
x = x0;
for k = 1:iters
	for i = 1:n
		x(i) = (1/A(i, i))*(b(i) - A(i, 1:n)*x + A(i, i)*x(i));
	end
end


