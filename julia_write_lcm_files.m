function julia_write_lcm_files(tmp_dir)
% Write the LCModel input files from a Julia reconstruction.
%
% run_julia_reco.jl reconstructs the CSI and writes CombinedCSI.mat next to the
% other results. This is the MATLAB side of that route: it assembles the arguments
% Write_LCM_files needs and produces the same .RAW, .control and batch files that
% MRSI_Reconstruction writes for the MATLAB route.
%
% The metadata is taken from Parameters.mat rather than from the Julia output, so
% both routes describe the spectra with the same header values. The dwell time is
% the exception: the Julia reconstruction's value is the one the spectra were sampled
% with, which MRSI.jl takes from the ADC dwell where the protocol states another.

load([tmp_dir '/Parameters.mat'],'Par')

if(~isfield(Par, 'ServerInfo') || ~isfield(Par.ServerInfo,'RunLCModel_CPUCores'))
    Par.ServerInfo.RunLCModel_CPUCores = 20;
end

[Data.csi, Dwelltime] = LoadJuliaCSI(Par);
Data = AddWaterReference(Data, Par, Dwelltime);
mask = LoadMask(tmp_dir, Par);

Paths.out_dir = [Par.Paths.out_path '/spectra'];
Paths.basis_file = Par.Paths.basis_path;
Paths.LCM_ProgramPath = Par.Paths.LCM_Path;
Paths.batchdir = tmp_dir;
if(Par.Flags.priors_flag)
    Paths.priors_dir = PriorsDirectory(Par);
end

MetaInfo.DimNames = {'x','y','z'};
MetaInfo.pat_name = Par.CSI.PatName;
MetaInfo.LarmorFreq = Par.CSI.LarmorFreq;
MetaInfo.dwelltime = Dwelltime;

if(isfield(Par.Paths,'LCM_ControlPath'))
    ControlInfo = Par.Paths.LCM_ControlPath;
else
    ControlInfo = 0;
end

Write_LCM_files(Data,Paths,MetaInfo,ControlInfo,mask,Par.ServerInfo.RunLCModel_CPUCores)

end



function [csi, Dwelltime] = LoadJuliaCSI(Par)
% run_julia_reco.jl stores the array as csi.Data, the field MRSI_Reconstruction
% uses, so the same variable serves both routes, and its dwell time in ns as
% csi.RecoPar.Dwelltimes.

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

Dwelltime = Par.CSI.Dwelltimes(1);
if(isfield(Loaded.csi,'RecoPar') && isfield(Loaded.csi.RecoPar,'Dwelltimes'))
    Dwelltime = double(Loaded.csi.RecoPar.Dwelltimes(1));
    if(abs(Dwelltime - Par.CSI.Dwelltimes(1)) > max(0.02 * Dwelltime, 5000))
        fprintf(['\nThe header reader gives a dwell time of %g ns, the Julia reconstruction ' ...
                 '%g ns. The LCModel files use the latter.\n'], Par.CSI.Dwelltimes(1), Dwelltime);
    end
end

end



function PriorsDir = PriorsDirectory(Par)
% The three prior maps -P points at, as raw files. Write_LCM_files reports a missing
% one through a file handle it never opens, so the check is made here.

PriorsDir = Par.Paths.priors_path;
for Name = {'0_pha_map', '1_pha_map', 'shift_map'}
    if(~exist([PriorsDir '/' Name{1} '.raw'],'file'))
        error('julia_write_lcm_files:MissingPrior', ...
              'The prior knowledge directory %s has no %s.raw.', PriorsDir, Name{1});
    end
end

end



function Data = AddWaterReference(Data, Par, Dwelltime)
% W1 hands LCModel the water FID of every voxel, which it scales the fit by.
% Write_LCM_files takes it as the array Data.watref on the grid of the spectra
% and zero-fills it to their length; the water pass left it in WaterReference.mat
% as WaterReferenceCSI, with Data and, from the Julia pass, RecoPar.Dwelltimes.

if(~(Par.Flags.WaterReference_flag && strcmp(Par.Settings.WaterReference_method,'W1')))
    return
end
WaterPath = [Par.Paths.out_path '/WaterReference.mat'];
if(~exist(WaterPath,'file'))
    error('julia_write_lcm_files:MissingWater', ...
          'W1 water referencing was asked for, but the water pass left no %s.', WaterPath);
end
% The MATLAB reconstruction leaves the water it gave LCModel, already on the grid and in
% the amplitude scale of the spectra. Write the same one again.
Loaded = load(WaterPath);
if(isfield(Loaded,'WaterReferenceOnCSIGrid'))
    if(~isequal(size(Loaded.WaterReferenceOnCSIGrid,1:3), size(Data.csi,1:3)))
        error('julia_write_lcm_files:WaterGrid', ...
              'WaterReferenceOnCSIGrid in %s is on a %s grid and the spectra on %s.', ...
              WaterPath, mat2str(size(Loaded.WaterReferenceOnCSIGrid,1:3)), mat2str(size(Data.csi,1:3)));
    end
    Data.watref = Loaded.WaterReferenceOnCSIGrid;
    return
end
if(~isfield(Loaded,'WaterReferenceCSI'))
    error('julia_write_lcm_files:BadWater', '%s carries no WaterReferenceCSI.', WaterPath);
end
Water = Loaded.WaterReferenceCSI;
if(~isfield(Water,'Data'))
    error('julia_write_lcm_files:BadWater', '%s carries no WaterReferenceCSI.Data.', WaterPath);
end

% LCModel reads the water file with the parameters of the spectra: same voxel, same
% dwell time. A water reference acquired on another grid or bandwidth has to be
% brought onto them first, which nothing here does.
WaterGrid = [size(Water.Data,1) size(Water.Data,2) size(Water.Data,3)];
CsiGrid = [size(Data.csi,1) size(Data.csi,2) size(Data.csi,3)];
if(~isequal(WaterGrid, CsiGrid))
    error('julia_write_lcm_files:WaterGrid', ...
          ['The water reference is on a %s grid and the spectra on %s. W1 scales voxel ' ...
           'by voxel, so the two have to share the grid.'], mat2str(WaterGrid), mat2str(CsiGrid));
end
% Write_LCM_files zero-fills a shorter water FID to the spectra's length but does not
% cut a longer one, and LCModel reads the file with the spectra's point count.
if(size(Water.Data,4) > size(Data.csi,4))
    error('julia_write_lcm_files:WaterLength', ...
          ['The water reference has %d points per voxel and the spectra %d. LCModel reads ' ...
           'the water file with the latter, so the longer one has to be cut first.'], ...
          size(Water.Data,4), size(Data.csi,4));
end
WaterDwell = [];
if(isfield(Water,'RecoPar') && isfield(Water.RecoPar,'Dwelltimes'))
    WaterDwell = Water.RecoPar.Dwelltimes(1);
elseif(isfield(Water,'Par') && isfield(Water.Par,'Dwelltimes'))
    WaterDwell = Water.Par.Dwelltimes(1);
end
% The two values come from different readers, one of which rounds to 10 us, so only
% a different bandwidth counts, not two roundings of the same one.
if(~isempty(WaterDwell) && abs(WaterDwell - Dwelltime) > max(0.02 * Dwelltime, 5000))
    error('julia_write_lcm_files:WaterDwell', ...
          ['The water reference was sampled with a dwell time of %g ns and the spectra ' ...
           'with %g ns. LCModel reads both with the latter, so the water area would be off.'], ...
          WaterDwell, Dwelltime);
end
Data.watref = Water.Data;

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
