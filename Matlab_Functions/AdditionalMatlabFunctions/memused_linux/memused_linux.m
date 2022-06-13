function [MemUsedByUser,MemFree,MemUsedTotal] = memused(quiet_flag)
%
% MemUsedByUser_linux Show memory usage of MATLAB and free memory in GiB.
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function shows the memory that is used by MATLAB in percent of all memory and the free system memory. Should work on all linux systems, probably even on all UNIX
% 
%
%
% MemUsedByUser = MemUsedByUser_linux(quiet_flag)
%
% Input: 
% -         quiet_flag                  ...    If 1, nothing is printed to display.
%
% Output:
% -         MemUsedByUser                     ...    used memory by MATLAB, in percent of all memory
% -         memfree                     ...    free memory on system
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions


% 0.1 Preparations

if(~exist('quiet_flag','var'))
    quiet_flag = 0;
end

% Check if we are on a linux system
[UnixSystem,MemFree] = unix('free -g');
UnixSystem = ~UnixSystem;
if(~UnixSystem)
	if(~quiet_flag)
		fprintf('\nWarning: memused cannot run on non-linux systems')
	end
	MemUsedByUser=0;MemUsedTotal = 0; MemFree=0;
	return;
end

% 0.3 Definitions
    






%% 1. Find out used memory


% Find out username.
[stat,uname] = unix('whoami');
[stat,uid] = unix('id');
uid = regexp(regexp(uid,'uid=\d+','match'),'\d+','match');
uid = uid{:}; uid = str2num(uid{:});

[stat,MemUsedByUser] = unix(['ps aux | grep "' uname '\|' uid '" | awk ''{print $4"\t"$11}'' | grep -i "matlab" | cut -f 1']);     % unix just performs the unix-command.
MemUsedByUser = sum(str2num(MemUsedByUser));                                                               % str2double does not work

[stat,MemUsedTotal] = unix('free -m | awk ''NR==3'' | awk ''{print $4}''');                         % unix just performs the unix-command. perform free -g command, take 3rd line of that, and 4th column.


[stat,MemFree] = unix('free -m | awk ''NR==3'' | awk ''{print $3}''');                         % unix just performs the unix-command. perform free -g command, take 3rd line of that, and 4th column.
clear stat
MemFree = sum(str2num(MemFree));                                                               % str2double does not work


%% 2. Display used memory

if(~quiet_flag)
    display([ char(10) num2str(MemUsedByUser) '% of total memory is used by MATLAB.' char(10)])
    display([ char(10) num2str(MemFree) ' MB of memory available.' char(10)])    
end




%% 3. Postparations

% fclose(fid)






