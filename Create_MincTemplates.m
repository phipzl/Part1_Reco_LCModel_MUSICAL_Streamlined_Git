%% -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                 PROGRAM TO HELL AND BACK                           %%%%%%%%%%%%%%
%%%%%%%%%%%%%%     ENTSTANDEN DURCH DIE KUNST DES PROGRAMMIERENS DURCH KONSEQUENTES ANSTARREN     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Create a figure template of a duerer image in the necessary size for rawtominc   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. DEFINITIONS, PREPARATIONS




%% 1. Create an image of a duerer image

duerer_orig = load('durer');
duerer_orig = duerer_orig.X;



%% 3. Resample it to the necessary size


zff = Par.Settings.ZeroFillMetMaps;
duerer = imresize(duerer_orig,[Par.CSI.nPhasEnc,Par.CSI.nFreqEnc],'bicubic');
duerer = repmat(duerer,[1 1 Par.CSI.nPartEnc*Par.CSI.nSLC]);
duerer(:,:,1:2:end) = flip(duerer(:,:,1:2:end),2);                % Flip every second slice, so that we can distinguish between the slices in the mincfile

duerer_zf = imresize(duerer_orig,[zff*Par.CSI.nPhasEnc,zff*Par.CSI.nFreqEnc], 'bicubic');
if(Par.CSI.ThreeD_flag)
     duerer_zf = repmat(duerer_zf,[1 1 zff*Par.CSI.nPartEnc*Par.CSI.nSLC]);
else
    duerer_zf = repmat(duerer_zf,[1 1 Par.CSI.nPartEnc*Par.CSI.nSLC]);   
end
duerer_zf(:,:,1:2:end) = flip(duerer_zf(:,:,1:2:end),2);
    


%% 4. WRITE THE FILES TO A .raw FILE


csi_template_fid = fopen([tmp_dir '/csi_template.raw'],'w');
csi_template_zf_fid = fopen([tmp_dir '/csi_template_zf.raw'],'w');

fwrite(csi_template_fid,duerer,'float');
fwrite(csi_template_zf_fid,duerer_zf,'float');

fclose(csi_template_fid);
fclose(csi_template_zf_fid);



%% 5 . Create image - Template

if(~Par.Flags.T1w_flag && isfield(Par,'Image'))
    duerer_image = imresize(duerer_orig,[Par.Image.nFreqEnc,Par.Image.nPhasEnc], 'bicubic');
    duerer_image = repmat(duerer_image, [1 1 Par.Image.nPartEnc*Par.Image.nSLC]);
    image_template_fid = fopen([tmp_dir '/image_template.raw'],'w');
	fwrite(image_template_fid,duerer_image,'float');       
    fclose(image_template_fid);
end



%% Orig (before user-set interpolation)

if(Par.Flags.InterpolateCSIResolution_flag == 1)
    duerer_BefInterpol = imresize(duerer_orig,[Par.CSI.nPhasEnc_BefInterpol,Par.CSI.nFreqEnc_BefInterpol], 'bicubic');
    duerer_BefInterpol = repmat(duerer_BefInterpol, [1 1 Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol]);
    Orig_template_fid = fopen([tmp_dir '/csi_template_BefInterpol.raw'],'w');
	fwrite(Orig_template_fid,duerer_BefInterpol,'float');       
    fclose(Orig_template_fid);
	
	bashstring_BefInterpol = ['rawtominc -clobber -float ' tmp_dir '/csi_template_BefInterpol.mnc -input ' tmp_dir '/csi_template_BefInterpol.raw'];
	bashstring_BefInterpol = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring_BefInterpol,Par.CSI.StepRead_BefInterpol, Par.CSI.StepPhase_BefInterpol, Par.CSI.StepSlice_BefInterpol);
	bashstring_BefInterpol = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring_BefInterpol, min(Par.CSI.POS_X_FirstVoxel_BefInterpol), ...
						min(Par.CSI.POS_Y_FirstVoxel_BefInterpol), min(Par.CSI.POS_Z_FirstVoxel_BefInterpol));
	bashstring_BefInterpol = sprintf('%s',bashstring_BefInterpol);
	bashstring_BefInterpol = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f', bashstring_BefInterpol, Par.CSI.ReadNormalVector, ...
	Par.CSI.PhaseNormalVector, Par.CSI.SliceNormalVector );            
	bashstring_BefInterpol = sprintf('%s %d %d %d',bashstring_BefInterpol,Par.CSI.nPartEnc_BefInterpol*Par.CSI.nSLC_BefInterpol,Par.CSI.nPhasEnc_BefInterpol,Par.CSI.nFreqEnc_BefInterpol);	
else
	bashstring_BefInterpol = '';
end


%% 6. Create Minc Files


% Normal
bashstring = ['rawtominc -clobber -float ' tmp_dir '/csi_template.mnc -input ' tmp_dir '/csi_template.raw'];
bashstring = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring,Par.CSI.StepRead, Par.CSI.StepPhase, Par.CSI.StepSlice);
bashstring = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring, min(Par.CSI.POS_X_FirstVoxel), min(Par.CSI.POS_Y_FirstVoxel), min(Par.CSI.POS_Z_FirstVoxel));
bashstring = sprintf('%s',bashstring);
bashstring = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f', bashstring, Par.CSI.ReadNormalVector, ...
Par.CSI.PhaseNormalVector, Par.CSI.SliceNormalVector );            
bashstring = sprintf('%s %d %d %d',bashstring,Par.CSI.nPartEnc*Par.CSI.nSLC,Par.CSI.nPhasEnc,Par.CSI.nFreqEnc);


% ZeroFilled
Pos_zf = [min(Par.CSI.POS_X_FirstVoxel)-Par.CSI.StepRead/2 + Par.CSI.StepRead_zf/2, ...
    min(Par.CSI.POS_Y_FirstVoxel)-Par.CSI.StepPhase/2+Par.CSI.StepPhase_zf/2, ...
    min(Par.CSI.POS_Z_FirstVoxel)-Par.CSI.StepSlice/2+Par.CSI.StepSlice_zf/2];
bashstring_zf = ['rawtominc -clobber -float ' tmp_dir '/csi_template_zf.mnc -input ' tmp_dir '/csi_template_zf.raw'];
bashstring_zf = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring_zf,Par.CSI.StepRead_zf, Par.CSI.StepPhase_zf, Par.CSI.StepSlice_zf);
bashstring_zf = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring_zf, Pos_zf(1), Pos_zf(2), Pos_zf(3)); % min correct?    
bashstring_zf = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f',bashstring_zf,Par.CSI.ReadNormalVector, ...
Par.CSI.PhaseNormalVector,Par.CSI.SliceNormalVector);    
if(Par.CSI.ThreeD_flag)
    bashstring_zf = sprintf('%s %d %d %d',bashstring_zf,zff*Par.CSI.nPartEnc*Par.CSI.nSLC,zff*Par.CSI.nPhasEnc,zff*Par.CSI.nFreqEnc);  
else                                                                                                                             
    bashstring_zf = sprintf('%s %d %d %d',bashstring_zf,Par.CSI.nPartEnc*Par.CSI.nSLC,zff*Par.CSI.nPhasEnc,zff*Par.CSI.nFreqEnc);  
end


bashstring_AlignFreq = '';
if(Par.Flags.AlignFreq_flag)
    if(isfield(Par.Paths,'AlignFreq_path') && isempty(regexpi(Par.Settings.AlignFreq_method,'Align')))
        duerer_AlignFreq = imresize(duerer_orig,[Par.AlignFreq.nPhasEnc,Par.AlignFreq.nFreqEnc], 'bicubic');
        duerer_AlignFreq = repmat(duerer_AlignFreq, [1 1 Par.AlignFreq.nPartEnc]);
        image_template_fid = fopen([tmp_dir '/AlignFreq_template.raw'],'w');
        fwrite(image_template_fid,duerer_AlignFreq,'float');       
        fclose(image_template_fid);

        bashstring_AlignFreq = ['rawtominc -clobber -float ' tmp_dir '/AlignFreq_template.mnc -input ' tmp_dir '/AlignFreq_template.raw'];
        bashstring_AlignFreq = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring_AlignFreq,Par.AlignFreq.StepRead, Par.AlignFreq.StepPhase, Par.AlignFreq.StepSlice);
        bashstring_AlignFreq = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring_AlignFreq, min(Par.AlignFreq.POS_X_FirstVoxel), min(Par.AlignFreq.POS_Y_FirstVoxel), min(Par.AlignFreq.POS_Z_FirstVoxel));
        bashstring_AlignFreq = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f', bashstring_AlignFreq, Par.CSI.ReadNormalVector, ...
        Par.CSI.PhaseNormalVector, Par.CSI.SliceNormalVector );            
        bashstring_AlignFreq = sprintf('%s %d %d %d',bashstring_AlignFreq,Par.AlignFreq.nPartEnc,Par.AlignFreq.nPhasEnc,Par.AlignFreq.nFreqEnc);
    else
        bashstring_AlignFreq = ['cp ' tmp_dir '/csi_template.mnc ' tmp_dir '/AlignFreq_template.mnc'];
    end
end



% Magnitude Template
if(~Par.Flags.T1w_flag)
    if(isfield(Par,'Image'))
        ImageOrCSI = 'Image';
        InterpolAddSuffix = '';
        Dimsize_image(3) = Par.Image.nPartEnc*Par.Image.nSLC;    % If Hadamard, but 1 slice-image for each CSI slice       
    else
        ImageOrCSI = 'CSI';
        if(isfield(Par.CSI,'StepRead_BefInterpol'))
            InterpolAddSuffix = '_BefInterpol';                                     % Only in the CSI case we need '_BefInterpol'
        else
            InterpolAddSuffix = '';            
        end
        Dimsize_image(3) = Par.CSI.nPartEnc*Par.CSI.nSLC;         
    end
    
    Step_image(1) = eval(['Par.' ImageOrCSI '.StepRead' InterpolAddSuffix]);
    Step_image(2) = eval(['Par.' ImageOrCSI '.StepPhase' InterpolAddSuffix]);
    Step_image(3) = eval(['Par.' ImageOrCSI '.StepSlice' InterpolAddSuffix]);
    Dimsize_image(1) = eval(['Par.' ImageOrCSI '.nFreqEnc' InterpolAddSuffix]);    
    Dimsize_image(2) = eval(['Par.' ImageOrCSI '.nPhasEnc' InterpolAddSuffix]);       
    
    bashstring_image = ['rawtominc -clobber -float ' tmp_dir '/mag_template.mnc -input ' tmp_dir '/' lower(ImageOrCSI) '_template' InterpolAddSuffix '.raw'];
    bashstring_image = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring_image, Step_image(1), Step_image(2),Step_image(3));    
    bashstring_image = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring_image, min(Par.CSI.POS_X_FirstVoxel), min(Par.CSI.POS_Y_FirstVoxel), min(Par.CSI.POS_Z_FirstVoxel));        
    bashstring_image = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f',bashstring_image,Par.CSI.ReadNormalVector,Par.CSI.PhaseNormalVector,Par.CSI.SliceNormalVector );            
    bashstring_image = sprintf('%s %d %d %d',bashstring_image,Dimsize_image(3),Dimsize_image(2),Dimsize_image(1)); 
    
else
    bashstring_image = '';
end

CopyCommands = sprintf('cp %s/csi_template.mnc %s/maps/csi_template.mnc', tmp_dir,Par.Paths.out_path);
CopyCommands = sprintf('%s\ncp %s/csi_template_zf.mnc %s/maps/csi_template_zf.mnc', CopyCommands, tmp_dir, Par.Paths.out_path);

fiddy = fopen([tmp_dir '/CreateMincTemplates.sh'],'w+');
fprintf(fiddy, '%s\n%s\n%s\n%s\n%s\n%s', bashstring, bashstring_zf,bashstring_BefInterpol, bashstring_image,bashstring_AlignFreq, CopyCommands);
fclose(fiddy);
fileattrib([tmp_dir '/CreateMincTemplates.sh'],'+x','gu');   % Executable for group and user.


% Why can't I run the commands from here in MATLAB ???????
% [bla,bla2] =  unix(bashstring);
% unix(bashstring_zf);
%
% Or at least running the .sh file
% ! bash ' tmp_dir '/CreateMincTemplates.sh


