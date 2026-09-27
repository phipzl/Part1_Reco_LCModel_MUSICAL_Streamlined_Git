function WaterData = water_reference_on_csi_grid(Water, csi)
% The W1 water FIDs on the grid of the metabolite CSI, for Write_LCM_files' Data.watref.
%
% Water is WaterReferenceCSI, the water pass's combined CSI; csi the metabolite CSI after
% its coil combination with the water pass's weights. LCModel scales every voxel by the
% water FID it is given, so the water has to describe the same voxel, with the same dwell
% time and not more points than the spectra.
%
% Amplitude: the shared coil weights cancel the reference scans, which leaves the two
% reconstruction gains. op_CalcAndApplyDensComp scales the density compensation of
% every partition so that a uniform image of the FoV reconstructs to a mean of 0.15, and
% the partitions are summed by an unnormalised FFT. Per unit spin density the gain is
% then 0.15 x FoV volume x the centre-to-mean ratio of that reconstructed uniform image,
% which the number of equidistant circles and their points leave unchanged to 0.1%. So
% scans that share the FoV share the gain, and the k-space resize below keeps voxel
% values.
%
% The resize assumes both scans share the FoV, position and orientation, which is
% checked. Other geometries need a resampling in world coordinates.

WaterPar = Water.RecoPar;
CsiPar = csi.RecoPar;

Fields = {'FoV_Read','FoV_Phase','FoV_Partition','Pos_Sag','Pos_Cor','Pos_Tra'};
for ii = 1:numel(Fields)
    if(abs(WaterPar.(Fields{ii}) - CsiPar.(Fields{ii})) > 0.1)       % mm
        error('water_reference_on_csi_grid:Geometry', ...
              ['The water reference has %s = %g mm and the spectra %g mm. W1 resizes the water in ' ...
               'k-space, which needs the same FoV, position and orientation.'], ...
              Fields{ii}, WaterPar.(Fields{ii}), CsiPar.(Fields{ii}));
    end
end
Normal = @(P) [P.SliceNormalVector_x P.SliceNormalVector_y P.SliceNormalVector_z];
if(norm(Normal(WaterPar) - Normal(CsiPar)) > 1e-3 || abs(WaterPar.InPlaneRotation - CsiPar.InPlaneRotation) > 1e-3)
    error('water_reference_on_csi_grid:Geometry', ...
          'The water reference and the spectra differ in orientation. W1 resizes the water in k-space, which needs the same one.');
end

% A dwell time rounded to 10 us by one reader and not by the other counts as the same.
if(abs(WaterPar.Dwelltimes(1) - CsiPar.Dwelltimes(1)) > max(0.02*CsiPar.Dwelltimes(1), 5000))
    error('water_reference_on_csi_grid:Dwell', ...
          'The water reference was sampled with a dwell time of %g ns and the spectra with %g ns.', ...
          WaterPar.Dwelltimes(1), CsiPar.Dwelltimes(1));
end
if(size(Water.Data,4) > size(csi.Data,4))
    error('water_reference_on_csi_grid:Length', ...
          'The water reference has %d points per voxel and the spectra %d. LCModel reads the water with the latter.', ...
          size(Water.Data,4), size(csi.Data,4));
end

% Same zero-filling or cutting in k-space that op_CoilCombineData applies to the weights.
CsiGrid = [size(csi.Data,1) size(csi.Data,2) size(csi.Data,3)];
WaterData = ZerofillOrCutkSpace(Water.Data, [CsiGrid size(Water.Data,4)], 1);

ReconstructionGainRatio = FoVVolume(CsiPar) / FoVVolume(WaterPar);
WaterData = WaterData * ReconstructionGainRatio;

fprintf('\nW1 water reference: %s resized to %s, reconstruction gain ratio %g.\n', ...
        mat2str(size(Water.Data,1:3)), mat2str(CsiGrid), ReconstructionGainRatio);

end



function V = FoVVolume(Par)
V = Par.FoV_Read * Par.FoV_Phase * Par.FoV_Partition;
end
