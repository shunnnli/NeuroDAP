function outputPath = osPathSwitch(input,options)

arguments
    input % string or char
    options.inputOS string = 'auto'
end

switchOS = false; 
server_mac = '/Volumes/MICROSCOPE/';
server_win = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\';

% Determine operation system
if ispc; targetOS = 'win';
elseif isunix; targetOS = 'mac';
end

% Detect input OS if needed
if strcmp(options.inputOS,'auto')
    if contains(input,'\') || contains(input,'research.files.med.harvard.edu')
        options.inputOS = 'win';
    elseif contains(input,'/') || contains(input,'Volumes')
        options.inputOS = 'mac';
    end
end

% Determine whether to swtich or not
if strcmp(options.inputOS,targetOS)
    outputPath = input;
    return
else
    switchOS = true;
end

% Detect and replace server path
if contains(input,server_mac) && switchOS % change from mac to windows
    splitted = strsplit(input,server_mac);
    outputPath = [server_win,splitted{2}];
    outputPath(strfind(outputPath,'/')) = '\';
elseif contains(input,server_win) && switchOS % change from windows to mac
    input(strfind(input,'\')) = '/';
    server_win(strfind(server_win,'\')) = '/';
    splitted = strsplit(input,server_win);
    outputPath = [server_mac,splitted{2}];
end

end