function outputPath = osPathSwitch(input,options)

arguments
    input % string or char
    options.inputOS string = 'auto'
end

% Determine operation system
if ispc; targetOS = 'win';
elseif isunix; targetOS = 'mac';
end

% Determine root path
if isunix
    if isfolder('/Volumes/Neurobio/MICROSCOPE'); rootPath = '/Volumes/Neurobio/MICROSCOPE';
    elseif isfolder('/Volumes/MICROSCOPE'); rootPath = '/Volumes/MICROSCOPE';
    else
        error('Did not contain root folder!');
    end
elseif ispc
    rootPath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE';
end

% Update root path
dirsplit = split(input,'MICROSCOPE');
outputPath = regexprep([rootPath,dirsplit{2}], '[\\/]', filesep);

end