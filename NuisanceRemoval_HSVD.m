% Remove nuisance (fat, water) signal from spectroscopy data using HSVD
% The functions used for this script are written by Chao Ma.


%% Some definitions


dim    = size(csi);
dt     = Par.CSI.Dwelltimes(1)*1e-9;
Factor = Par.CSI.LarmorFreq * 1e-6;

% run NuisanceRemovalControlPath. Please don't write crap in there that causes MATLAB to crash or worse...
run(Par.Paths.NuisRem_ControlPath)



% Signal selective HSVD params
selParams               = struct('dt',dt,'n',NoOfSingVals);                           % n: Number of singular values that should be kept in the truncated SVD.
optWat                  = struct('fSel2',(4.65-WaterPPMs)*Factor,'maxT2',WaterT2s); % Search for two signal ranges: One between -500 and 100 Hz and with a T2<1e6, and another between 100 and 120 Hz with a T2<20
optLip                  = struct('fSel2',(4.65-LipidPPMs)*Factor,'maxT2',LipidT2s);
optMeta.fSel2           = (4.65-MetaboPPMs)*Factor;                                        % Parameters of the metabolites: Any T2 btw freq 110 and 300 Hz. What is this used for?!
optOther                = [];%struct('fSel2',[-500, 500],'maxT2',[1e6]);    % If voxel is not in mask, process data with that settings. If empty, skip voxel.
selParams.signalType    = 'nuisance';
selParams.NtEcho        = 0;

%% water removal 


fprintf('\n\nRemove nuisance signals . . . ');
NuisTic = tic;
for i = 1:NuisRemIterations         % 3 for Phantom, 1 invivo
    for slc = 1:size(csi,3)
        [rho_W(:,:,1,:),csi(:,:,slc,:)] = nsRm(squeeze(csi(:,:,slc,:)), mask(:,:,slc), zeros(dim(1:2)), selParams,...
                                          optWat, optLip, optMeta, optOther);
    end
end
fprintf('%s seconds.',num2str(toc(NuisTic)))




