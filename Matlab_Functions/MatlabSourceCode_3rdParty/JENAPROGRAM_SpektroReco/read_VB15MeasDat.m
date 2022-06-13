function [data, prot, mdh] = read_VB15MeasDat(measDat)
%
% function [data, prot] = read_iceMeasDat(measDat);
%
%   This function read ICE meas.dat (VB13) files
%   Output is an n-dimensional complex array with the raw data
%   and a struct with the protocol parameters
%

% *************************************************************************

if nargin < 1                       % check for imput of meas file name 
    cur_dir =pwd;
    [n,pf]  = uigetfile('*.dat','Select your MeasOut-File!');
    if length(n) < 2 && length(pf) < 2
        data    = 0;
        prot    = 0;
        cd(cur_dir);
        return;
    end
    measDat     = strcat(pf,n);
    cd (cur_dir);
end

for i = length(measDat):-1:1
    if strcmp(measDat(i),'/') || strcmp(measDat(i),'\')
        name_meas = measDat(i+1:length(measDat));
        break;
    end
end

% *************************************************************************
% open file to read
fid_dat = fopen(measDat,'r');
ascLen  = fread(fid_dat, 1, 'int32');
%ascLen = ascLen/4;

% generate asc file and open it
ascFile = [measDat(1:(length(measDat)-3)) 'asc'];
fid_asc = fopen(ascFile, 'wb');

protdata = fread(fid_dat, ascLen, 'int8');
n = fwrite(fid_asc, protdata, 'int8');

fclose(fid_dat);
fclose(fid_asc);

% *************************************************************************
[data, prot, mdh]       = read_iceMeasOut(measDat, ascFile);