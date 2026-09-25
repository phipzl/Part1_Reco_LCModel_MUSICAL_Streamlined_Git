function julia_write_lcm_files(tmp_dir)
% Write the LCModel input files from a Julia reconstruction.
%
% run_julia_reco.jl reconstructs the CSI and writes CombinedCSI.mat next to the
% other results. This is the MATLAB side of that route: it assembles the arguments
% Write_LCM_files needs and produces the same .RAW, .control and batch files that
% MRSI_Reconstruction writes for the MATLAB route.
%
% The metadata is taken from Parameters.mat rather than from the Julia output, so
% both routes describe the spectra with the same header values.

load([tmp_dir '/Parameters.mat'],'Par')

if(~isfield(Par, 'ServerInfo') || ~isfield(Par.ServerInfo,'RunLCModel_CPUCores'))
    Par.ServerInfo.RunLCModel_CPUCores = 20;
end

Data.csi = LoadJuliaCSI(Par);
mask = LoadMask(tmp_dir, Par);

Paths.out_dir = [Par.Paths.out_path '/spectra'];
Paths.basis_file = Par.Paths.basis_path;
Paths.LCM_ProgramPath = Par.Paths.LCM_Path;
Paths.batchdir = tmp_dir;

MetaInfo.DimNames = {'x','y','z'};
MetaInfo.pat_name = Par.CSI.PatName;
MetaInfo.LarmorFreq = Par.CSI.LarmorFreq;
MetaInfo.dwelltime = Par.CSI.Dwelltimes(1);

if(isfield(Par.Paths,'LCM_ControlPath'))
    ControlInfo = Par.Paths.LCM_ControlPath;
else
    ControlInfo = 0;
end

Write_LCM_files(Data,Paths,MetaInfo,ControlInfo,mask,Par.ServerInfo.RunLCModel_CPUCores)

end



function csi = LoadJuliaCSI(Par)
% run_julia_reco.jl stores the array as csi.Data, the field MRSI_Reconstruction
% uses, so the same variable serves both routes.

CombinedPath = [Par.Paths.out_path '/CombinedCSI.mat'];
if(~exist(CombinedPath,'file'))
    error('julia_write_lcm_files:MissingInput', ...
          ['The Julia reconstruction output %s does not exist. ' ...
           'Check whether run_julia_reco.jl completed.'], CombinedPath);
end

Loaded = load(CombinedPath,'csi');
if(~isfield(Loaded,'csi') || ~isfield(Loaded.csi,'Data'))
    error('julia_write_lcm_files:BadInput', ...
          '%s does not contain csi.Data.', CombinedPath);
end
csi = Loaded.csi.Data;

end



function mask = LoadMask(tmp_dir, Par)
% Same mask and same fallback as MRSI_Reconstruction, so both routes fit the
% same voxels.

MaskSize = [Par.CSI.nFreqEnc Par.CSI.nPhasEnc Par.CSI.nSLC*Par.CSI.nPartEnc];
MaskPath = [tmp_dir '/mask_brain.raw'];

if(~exist(MaskPath,'file'))
    mask = ones(MaskSize);
    return
end

fid_mask = fopen(MaskPath,'r');
mask = reshape(fread(fid_mask, 'float'), MaskSize);
fclose(fid_mask);

if(sum(mask(:)) == 0)
    mask = ones(MaskSize);
end

end
