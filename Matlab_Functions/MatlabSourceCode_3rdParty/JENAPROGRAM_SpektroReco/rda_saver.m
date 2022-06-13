function rda_saver(data_rda,header_rda,name_rda,modus,machine_format)

if      nargin < 5 && nargin > 3
    machine_format  = 'l';
elseif  nargin < 4 && nargin > 2
    modus           = 'double';
    machine_format  = 'l';
elseif  nargin < 3 && nargin > 1
    cur_dir = pwd;
    [n,pf]  = uiputfile('*.rda','Erstellen eines rda-Files');
    if length(n) < 2 && length(pf) < 2
        cd(cur_dir);
        return;
    end
    name_rda        = strcat(pf,n);
    cd(cur_dir);    
    modus           = 'double';
    machine_format  = 'l';
end


% *************************************************************************

if iscell(data_rda)
    [N_CSI_slices, N_CSI_rows, N_CSI_colums]    = size(data_rda);
    N_samp  = max(size(data_rda{1,1,1}));
    temp    = zeros(1,N_CSI_slices*N_CSI_rows*N_CSI_colums*N_samp*2);
    start   = 1;
    for i = 1:N_CSI_slices
        for j = 1:N_CSI_rows
            for k = 1:N_CSI_colums
                temp1                           = data_rda{i,j,k}(:);
                temp2                           = zeros(1,2*N_samp);
                temp2(1:2:2*N_samp)             = real(temp1);
                temp2(2:2:2*N_samp)             = imag(temp1);
                temp(start:start+2*N_samp-1)    = temp2;
                start                           = start+2*N_samp;
                clear temp1 temp2
            end
        end
    end
    data_rda    = temp; clear start
else
    if ~isreal(data_rda)
        temp(1:2:2*length(data_rda))  = real(data_rda);
        temp(2:2:2*length(data_rda))  = imag(data_rda);
        data_rda    = temp;
    end
end

% *************************************************************************
fid = fopen(name_rda,'w');

fseek(fid, 0, -1);
fwrite(fid,header_rda,'char');
fwrite(fid,data_rda,modus,machine_format);

fclose(fid);