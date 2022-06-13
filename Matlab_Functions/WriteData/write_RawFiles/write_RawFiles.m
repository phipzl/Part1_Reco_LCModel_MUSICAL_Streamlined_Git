function write_RawFiles(RawArray,RawFile,precision)
%
% write_RawFiles Write raw-files.
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function can really do nothing, and more specifically, exactly nothing.
% 
%
%
% write_RawFiles(RawArray,RawFile,precision)
%
% Input: 
% -         inputvar1                   ...    This is the first input
% -         inputvar2                   ...    And this the second
%
% Output:
% -         None.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 





%% 0. Preparations


if(~exist('precision','var'))
	precision = 'float32';
end

Directory = regexp(RawFile,'/'); Directory = RawFile(1:Directory(end)-1);
if(~exist(Directory,'dir'))
	mkdir(Directory);
end




%% 1. Write data

raw_fid = fopen(RawFile,'w');
for Slice_no = 1:size(RawArray,3)
    fwrite(raw_fid, RawArray(:,:,Slice_no), precision);
end
fclose(raw_fid);

