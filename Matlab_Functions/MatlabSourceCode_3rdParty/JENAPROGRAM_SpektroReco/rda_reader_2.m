function [d,h] = rda_reader_2(read_flag,rda_name)

% read_flag define what should be readed:
%   read_flag = 0       - read rda header only, 'd' is then a zero matrix
%   read_flag = 1       - read rda header and data

cur_dir     = pwd;
if nargin == 1
    cur_dir     = pwd;
    [n,pf]      = uigetfile('*.*','Select rda-file!');
    if length(n) < 2 && length(pf) < 2
        cd(cur_dir);
        d           = 0;
        h           = 0;
        return;
    end
    rda_name    = strcat(pf,n);
    cd (cur_dir);
elseif  nargin < 1
    read_flag   = 1;
    cur_dir     = pwd;
    [n,pf]      = uigetfile('*.*','Select rda-file!');
    if length(n) < 2 && length(pf) < 2
        cd(cur_dir);
        d           = 0;
        h           = 0;
        return;
    end
    rda_name    = strcat(pf,n);
    cd (cur_dir);
end

% *************************************************************************
% Read rda header only

fid         = fopen(rda_name,'r');

test        = fread(fid, 16,'*char' );
if isempty(strfind(test',char(zeros(1, 16))))
    fseek(fid,0,-1);
    flag_stop   = 0;
    k_fid       = 0;
    while ~flag_stop
        fseek(fid,k_fid,-1);
        if strcmp('>>> End of header <<<',fread(fid, max(size('>>> End of header <<<')),'*char' )')

            N_header            = ftell(fid);
            fseek(fid,0,-1);
            header              = fread(fid,N_header,'*char');
            h                   = read_header(header');
            start_data          = ftell(fid)+2;
            flag_stop           = 1;
        end
        k_fid   = k_fid+1;
    end

    % *************************************************************************
    % extract header infos
    N_header   = max(size(h));

    for k = 1:N_header
        if strcmp(h{k,1},'VectorSize'),         N_samp          = str2double(h{k,2});   end
        if strcmp(h{k,1},'CSIMatrixSize[0]'),   CSI_size(1)     = str2double(h{k,2});   end
        if strcmp(h{k,1},'CSIMatrixSize[1]'),   CSI_size(2)     = str2double(h{k,2});   end
        if strcmp(h{k,1},'CSIMatrixSize[2]'),   CSI_size(3)     = str2double(h{k,2});   end
        % Definition: CSI_size == [N_CSI_colums N_CSI_rows N_CSI_slices]
    end

    d       = cell(CSI_size(3), CSI_size(2), CSI_size(1));
    for k_slice = 1:CSI_size(3)
        for k_row = 1:CSI_size(2)
            for k_col = 1:CSI_size(1)
                d{k_slice,k_row,k_col}  = zeros(N_samp,1);
            end
        end
    end

    if read_flag == 0
        fclose(fid);
        return,
    end

    % *************************************************************************
    % extract header infos and data

    if read_flag == 1
        fseek(fid,start_data,-1);
        full_data   = fread(fid,2*N_samp*CSI_size(1)*CSI_size(2)*CSI_size(3),'double');

        if max(size(full_data)) < 2*N_samp*CSI_size(1)*CSI_size(2)*CSI_size(3)
            fseek(fid,start_data-2,-1);
            full_data   = fread(fid,2*N_samp*CSI_size(1)*CSI_size(2)*CSI_size(3),'double');
        end

        start       = 1;

        k1_slice    = 1;
        for k_slice = 1:CSI_size(3)
            k1_row  = 1;
            for k_row = 1:CSI_size(2)
                k1_col  = 1;        
                for k_col = 1:CSI_size(1)
                    temp            = full_data(start:start+2*N_samp-1);
                    d{k1_slice, k1_row,k_col}       = complex(temp(1:2:2*N_samp),temp(2:2:2*N_samp));
                    start           = start+2*N_samp;
                    k1_col          = k1_col+1;
                end
                k1_row  = k1_row+1;
            end
            k1_slice  = k1_slice+1;
        end

        fclose(fid);
    end
else
    fclose(fid);
    cd(cur_dir);
    d           = 0;
    h           = 0;
    return;
end