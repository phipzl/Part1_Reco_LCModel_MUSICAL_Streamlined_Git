function header = read_header(h)

% h = original rda-header
% header = encoded rda-header as a cell-array
% pos_info = cell-arrary, with information about position of spectroscopy-voxel

% *************************************************************************
% read and transcript rda-header
k       = 1;
start   = 1;

for i = 1:length(h)-2
    if strcmp(h(i:i+1),[char(13),char(10)])
        ende        = i-1;
        temp        = h(start:ende);
        for j = 1:length(temp)-1
            if strcmp(temp(j:j+1),': ')
                header{k,1} = temp(1:j-1);
                header{k,2} = temp(j+2:length(temp));
            end
        end        
        start       = i+2;
        k           = k+1;
    end
end