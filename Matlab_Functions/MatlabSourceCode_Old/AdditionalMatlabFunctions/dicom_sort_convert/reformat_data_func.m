function data = reformat_data_func(data)

%   if the 'Advanced' feature are selected in dicom_sort_convert_main, this
%   function reformats a dataset with 3 or 4 dimensions into one with up to
%   7.  The dimensions of these are specified by the dims parameters, and
%   are, in order
%   [x][y][z][time or diffusion gradient][echos][channels][phase/mag]
%   The first 3 dimensions are always image x,y,z. After that, dimensions
%   are only used if variation in that dimension is present (more than one
%   time point, echo, channel, p/m)


% NON-MOSIAC
% odims(1) = x
% odims(2) = y
% odims(3) = channels-slices-echoes-(p/m in VA or sets?)-repetitions
 
% MOSIAC
% odims(1) = x
% odims(2) = y
% odims(3) = z
% odims(4) = slices-(p/m in VA or sets?)-channels-echo-repetitions

% dims(1) = x
% dims(2) = y
% dims(3) = z
% dims(4) = time
% dims(5) = echos
% dims(6) = channels
% dims(7) = phase/magnitude

dims(1:7) = 1;
odims = 0;
reformat = 'no';
datalabel = '[x,y,z';

convert_format = data.convert_format;
current_dir = data.current_dir;
current_nfiles = data.current_nfiles;
n_channels = data.current_header_pars.n_channels;

if data.current_analyze_struct.hdr.dime.dim(1) == 4
    odims (1:4) = data.current_analyze_struct.hdr.dime.dim(2:5);
else
    odims (1:3) = data.current_analyze_struct.hdr.dime.dim(2:4);
end

z = data.current_header_pars.dim_nslices;

dims (1) = odims(1);
dims (2) = odims(2);
dims (3) = z;

%   dim4 is for time (for EPI-type) or diffusion gradient
if data.current_header_pars.dim_nr > 1
    dims(4) = data.current_header_pars.dim_nr;    
    switch data.current_header_pars.sequence_type
        case 'ep2d_diff'
            if data.current_header_pars.dim_ngd > 1
                dims(4) = data.current_header_pars.dim_ngd;
                datalabel = [datalabel ',diff_g'];
                %DTI with multiple repetitions
                if (data.current_header_pars.dim_nr > 1)
                    disp('!!! lRepetitions gives NR > 1 - assuming this is number of diffusion directions');
                    dims(4) = data.current_header_pars.dim_nr;
                    data.warning_present = 'yes';
                    this_warning_text = 'lRepetitions gives NR > 1 - assuming this is number of diffusion directions - check reformatted data';
                    data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
                end
            end
        otherwise
            %   for time-series data
            datalabel = [datalabel ',t'];
    end
end

%   dim5 is echo number: true multi-echo
if (data.current_header_pars.nechos > 1 && strcmp(data.true_multi_echo,'yes')==1)
    dims(5) = data.current_header_pars.nechos;
    datalabel = [datalabel ',echo'];    
end

%   dim6 is channel
if strcmp(data.current_header_pars.sep_channels, 'yes')
    dims(6) = n_channels;    
    datalabel = [datalabel ',channel'];
end

%   dim7 is phase/magnitude (VA..., in VB these are in different scans)
if strcmp(data.current_header_pars.data_type, 'MP')
    dims(7) = 2;    
    datalabel = [datalabel ',phase/mag'];
end
datalabel = [datalabel ']'];

%   reformat if more than 4D, or if 4D and it is not a time-series
if size(dims(dims(4:end)~=1),2) > 1 || (size(dims(dims(4:end)~=1),2)==1 && dims(4)==1)
    reformat = 'yes';
end

%   below, some scans that shouldn't be reformatted - reverse decision
try
    if numel(findstr(data.current_dcminfo.ImageType, 'T2 MAP'))==1 ||  numel(findstr(data.current_dcminfo.ImageType, 'T2_STAR MAP'))==1
        reformat = 'no';
    end
catch
    disp('!!! In reformat_data_func.m: data.current_dcminfo.ImageType is not defined. Not reformatting');
    reformat = 'no';
end


if strcmp(reformat,'yes')
    disp([' - scan ' data.current_scan ' is ' datalabel ', reformatting into those dimensions']);
    current_reform_dir = fullfile(current_dir,'reform');
    s = warning ('query', 'MATLAB:MKDIR:DirectoryExists') ;	    % get current state
    warning ('off', 'MATLAB:MKDIR:DirectoryExists') ;
    mkdir(current_reform_dir);
    warning (s) ;						    % restore state
    switch convert_format
        case 'nifti'
            extension = '.nii';
        case 'analyze'
            extension = '';
    end
    readfile = fullfile(current_dir, sprintf('Image%s', extension));
    nfiles_required = prod(odims(3:end));
    if strcmp(data.current_mosaic_flag, 'mosaic') == 1
        nfiles_required = prod(odims(3:end))/z;
    end
    %   check that the product of the reformatted dimensions equals the old dimensions, and the number of files
    if nfiles_required ~= current_nfiles
        disp(sprintf('!!! The scan dimensions are %i, which means there should be %i files, but there are %i', num2str(odims), nfiles_required, current_nfiles));
        data.warning_present = 'yes';
        this_warning_text = sprintf('!!! There should be %i files in this series, but there are %i',prod(odims(3:end)), current_nfiles);
        data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
        return;
    end
    if prod(odims(3:end)) ~= prod(dims(3:end))
        dim_ratio = prod(odims(3:end))/prod(dims(3:end));
        disp(sprintf('!!! The old scan dimensions are %s, and the new dimensions %s, which differ by a factor %f', num2str(odims(1:end)), num2str(dims(1:end)), dim_ratio));
        if round(dim_ratio) == dim_ratio
            disp(sprintf('!!! As this is an integer, trying to continue, putting scans into a further dimension: check images', n_channels/dim_ratio, n_channels));
            data.warning_present = 'yes';
            this_warning_text = sprintf('!!! Added a further dimension in reformatting, with size %i ', dim_ratio);
            data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
            new_dim = dim_ratio;
            dims(current_dim) = new_dim;
            current_dim = current_dim + 1;
        else
            disp(sprintf('!!! The old scan dimensions are %s, and the new dimensions %s, which differ by the a factor %f', num2str(odims), num2str(dims), dim_ratio));
            data.warning_present = 'yes';
            this_warning_text = sprintf('!!! The number of old and new dimensions doesn''t match (scan incomplete?), not reformatting', n_channels, n_channels/dim_ratio);
            data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
            return;
        end
    end
    try
        new_image = zeros(dims,data.current_header_pars.precision);
    catch
        disp('!!! Could not create the image for reformatting');
        return;
    end
    try
        old_image_nii = load_nii(readfile);
    catch
        disp('!!! Could not load the old image to perform reformatting');
        data.warning_present = 'yes';
        this_warning_text = '!!! Image couldn''t be loaded prior to reformatting';
        data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
        return;
    end
    try
        switch data.current_mosaic_flag
            case 'mosaic'
                % The order of 'slices' in odims is (x, y) then slices, (sets?), channels, echoes, repetitions.
                new_image=reshape(old_image_nii.img, [dims(1) dims(2) dims(3) dims(7) dims(6) dims(5) dims(4)]);
                new_image=permute(new_image, [1 2 3 7 6 5 4]);
            otherwise
                % The order of 'slices' in odims is (x, y) channels, slices, echoes, (sets?), repetitions.
                new_image=reshape(old_image_nii.img, [dims(1) dims(2) dims(6) dims(3) dims(7) dims(5) dims(4)]);
                new_image=permute(new_image, [1 2 4 7 6 3 5]);
        end
        data.warning_present = 'yes';
        this_warning_text = ['was reformatted to ' datalabel ' in /reformat)' ];
        data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
        % sometime the slices are swapped in sep channel data
%         if strcmp(data.current_header_pars.sep_channels,'yes') && (data.current_header_pars.dim_nr > 1) && strcmp(data.current_mosaic_flag, 'mosaic') ~= 1
%             disp('  -  swapping slice order');
%             %%
%             if iseven(z)
%                 slice_acq_order=[2:2:z 1:2:z-1];
%             else
%                 slice_acq_order=[1:2:z 2:2:z];
%             end
%             [dummy, new_slice_order] = sort(slice_acq_order);
%             new_image=new_image(:,:,new_slice_order,:,:,:);
%             %%
%         end
    catch
        disp(['!!!' data.current_scan ': Couldn''t reformat data into new matrix - in reformat_data_func.m']);
        data.warning_present = 'yes';
        this_warning_text = ['Warning!!!: couldn''t be reformatted to ' datalabel];
        data.warning_text = sprintf('%sScan: %s. %s\n', data.warning_text, data.current_scan, this_warning_text);
end

old_image_hdr = load_nii_hdr(readfile);

new_image = squeeze(new_image);
size_new_image=size(new_image);
%   If there was only 1 slice, put this back in as singleton third dimension
if dims(3)==1
    new_image=reshape(new_image, [size_new_image(1:2) 1 size_new_image(3:end)]);
end
new_image_nii = make_nii_sr(new_image, 4);
new_image_nii.hdr = centre_header(new_image_nii.hdr);
new_image_nii.hdr.dime.pixdim(2:4) = old_image_hdr.dime.pixdim(2:4);
save_nii(new_image_nii, fullfile(current_reform_dir, ['Image' extension]));
end

