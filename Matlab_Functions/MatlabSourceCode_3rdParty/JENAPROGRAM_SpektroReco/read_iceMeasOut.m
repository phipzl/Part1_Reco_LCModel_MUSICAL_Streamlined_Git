function [data, prot, mdh] = read_iceMeasOut(measOut, measAsc)

% call the function to read asc file with the meat information corresponding to the current meas file
[prot]              = read_asc(measAsc);

% initialize the size of spectroscopic data matrix
nCha                = prot.lNumberOfChannels;       % Number of coil channels (is allways zero in prot structure, would be therefore adapted through the data reading)
nRep                = prot.lRepetitions+1;          % Number of measuremenst (for dynamic runs)
nAve                = prot.lAverages;               % Number of averages (NAS)
N_samp              = prot.svs_vector_size;         % Number of complex data points in FID (is allways a half of initial vector length in meas data) 

% create the zero data matrix
data                = zeros(nRep, nAve, nCha, 2*N_samp);

% *************************************************************************
% create file ID and open the meas.out file
fid                 = fopen(measOut, 'r');

% length of dummy header in front of meas data (depends on the NUMARIS version)
if regexpi(prot.NumarisVersion, 'N4_VB[0-9]+')              % N4_VB15 or 17 are included here because of [0-9]+
    dummyLength = fread(fid, 1, 'int32');
    fread(fid, dummyLength-4, 'int8');
else
    dummyLength = 32;
    fread(fid, dummyLength, 'int8');
end;

% estimate the file length
cur_fid         = ftell(fid);
fseek(fid,0,1);
end_fid         = ftell(fid);
fseek(fid,cur_fid,-1);

% *************************************************************************
% read the meas data

meta            = zeros(128,2);    % WHY 128?
k               = 1;

while ftell(fid) < end_fid
    mdh                                                 = read_mdh(fid);        % read the header (mdh) of the current data segment (corresponds to following variable 'tmp_line')
    [tmp_line,N]                                        = fread(fid, mdh.ushSamplesInScan*2, 'single');     % read the current data segment; single = float = float32; [A,count] = fread(...), count = number of elements read into A
    
    if N < 0.1*N_samp       % check whether the current data segment is data or end of file information 
        break;              % end the data reading, if the length of current data segment is smaller than the expected length of a spectroscopic data vector  
    else
        
        chan                                                = mdh.ulChannelId+1;        % coil channel ID
        set                                                 = mdh.sLC.ushSet+1;         % average ID
        repet_counter                                       = mdh.sLC.ushEcho+1;        % dynamic repetition ID or editing ID, e.g. in MEGAPRESS data 
        
        meta(k,1)       = set; meta(k,2)       = chan;  k   = k+1;
        
        % feel the data matrix        
        data(repet_counter,set,chan,1:length(tmp_line)/2)  = complex(   tmp_line(1:2:length(tmp_line)), ...     % real part of complex spectroscopic data
                                                                        tmp_line(2:2:length(tmp_line)));        % imaginary part of complex spectroscopic data
    end
end

fclose(fid);