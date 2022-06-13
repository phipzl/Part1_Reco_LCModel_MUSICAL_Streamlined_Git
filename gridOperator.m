
function [H,APF,DCF] = gridOperator(nSlices,slicenmb,Traj,M,overgrid,kernel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [H,ACF,DCF] = gridOperator(Traj,M,overgrid,kernel)
%
% INPUT: 2D Traj            (2 x Nsamples) defined between -0.5 and 0.5
%        M                  integer; ImgSize without oversampling
%        overgrid           overgridding by factor os
%        kernel             structure of the convolution kernel with
%
%                           kernel.values
%                           kernel.kwidth     (kwidth = width/2/M)
%
%
%
% OUTPUT: H                 Sparse Matrix; Gridding Operator   (Nsamples x N^2)
%         ACF               Apodization Compensation Function  (NxN)
%         DCF               Density Compensation Function      (Nsamples x 1)
%
%                           N = ceil(M*overgrid/2)*2  
%         
%         -----------------------------------------------------------------------------
%         S                 kernel values to build spares matrix
%         I                 index of non-Cartesian data point
%         J                 index of position in Cartesian K-space (J = (y-1)*N + x))
%
%         H = sparse(I,J,S,Nsamples,N*N,length(S));
%         -----------------------------------------------------------------------------
%
% written by Felix Breuer Fraunhofer IIS/MRB 23.02.2016
% Based on Code from Brian Hargreaves
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check Inputs

if nargin < 2
    error('At least 2 Inputs are required')
end

if mod(M,2)
    warning('Input M should be even')
end

%% Use default overgridding
if ~exist('overgrid','var');
    overgrid = 2;
else
    if (numel(overgrid) ~= 1 || overgrid > 5 || overgrid < 1),
        error('Check overgridding factor: should be 1 < os < 5'), end
end

%N = ceil(M*overgrid/2)*2;
N = M; %LUKAS: HACKERINO


%% Use default Kernel
if ~exist('kernel','var');
%     display('Using default Gridding Kernel (Kaiser Bessel)')
    width = 3;       % Gridding kernel width in dk without os
    ns    = 100;     % Should be greater than 32
    beta  = pi*sqrt(width^2/overgrid^2*(overgrid-0.5)^2-0.8);	% beta by Beatty et al.
    vals  = kaiser(2*ns,beta);  vals(1:ns) = [];
    % Setup the kernel structure for the convolution kernel:
    kernel.values = vals;
    kernel.kwidth = width/2/M;
else
    if ~all(isfield(kernel,{'values','kwidth'})),
        error('Check kernel structure');
    end
    
    if (kernel.kwidth*M > M/4);
        error('Check kernel width:');
    end
    
    if (any(kernel.values < 0) || numel(kernel.values(:)) < 32);
        error('Check kernel values: values should be positive and number of values should be greater than 32 ');
    end
    
end

%% Set up gridding operator

Traj = Traj(:,:);

nKernelpoints = length(kernel.values);
kwidth        = kernel.kwidth;
kerneltable   = kernel.values;
nSamples      = size(Traj(:,:),2);

% Calculate the kernel extend on Cartesian grid for each sampled
% non-Cartesian location

iMin = floor((Traj-kwidth)*N + N/2+1);
iMax = floor((Traj+kwidth)*N + N/2+1);

% Allocate memory for variables in sparse matrix

nKernelVals = round((kwidth*N).^2*pi);

I = zeros(nSamples*nKernelVals,1);
J = zeros(nSamples*nKernelVals,1);
S = zeros(nSamples*nKernelVals,1);

cnt = 0;

for k = 1:nSamples,
    
    for x = iMin(1,k):iMax(1,k)
        dkx = (x - (N/2+1))/N - Traj(1,k);
        
        for y = iMin(2,k):iMax(2,k)
            dky = (y - (N/2+1))/N - Traj(2,k);
            
            dk = sqrt(dkx.^2+dky.^2);
            
            if (dk<kwidth)
                cnt = cnt+1;
                %Find index in kernel lookup table
                ind =  floor(dk/kwidth*(nKernelpoints-1))+1;
                kern = kerneltable(ind);
                % circular boundary condition
                y_ = mod(y-1,N)+1;
                x_ = mod(x-1,N)+1;
                
                J(cnt) = x_ + N*(y_-1);
                I(cnt) = k;
                S(cnt) = kern;
                
                
            end
        end
        
    end
end

I(cnt+1:end) = [];
J(cnt+1:end) = [];
S(cnt+1:end) = [];

% Convolution Kernel in sparse Matrix form
H = sparse(I,J,S,nSamples,N*N,length(S));
% Normalization
H = H./sum(full(H(1,:)));


%% Appodization Function
if nargout > 1
    W = abs(ifftshift(ifft2(ifftshift(reshape(full(H(1,:)),[N,N])))));
    W = W./max(W(:));
    W(W<1e-2) = 1e-2;
    APF = 1./W;
end

HAM=0.54+0.46*cos(2*pi*sqrt(Traj(1,:).^2+Traj(2,:).^2));
%HAM=0.54+0.46*cos(2*pi*sqrt(Traj(1,:).^2+Traj(2,:).^2+((nSlices+1)/2-slicenmb)^2));
%HAM'.*


%% Density Compensation Function: Pipe et al.
if nargout > 2
    DCF = ones(nSamples,1);
    niter = 5;
    for i = 1:niter;
        goal = H*(H'*DCF);
        DCF = HAM'.*DCF./abs(goal);
    end
    DCF = reshape(DCF,[size(Traj,2),size(Traj,3)]);
end



end


