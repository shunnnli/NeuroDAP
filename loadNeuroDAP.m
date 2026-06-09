function neuroDAPDir = loadNeuroDAP()
% loadNeuroDAP  Add the NeuroDAP Methods folder to the MATLAB path.
%
%   neuroDAPDir = loadNeuroDAP()
%
%   This function automatically finds the main NeuroDAP folder by starting
%   from the folder of the script/function that called loadNeuroDAP, then
%   walking upward until it finds a folder that contains "Methods". If MATLAB
%   runs the caller from a temporary Editor_* folder, it also tries the folder
%   containing loadNeuroDAP.m and the current working folder.
%
%   This avoids hard-coded paths such as:
%       /Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods
%
%   It should work across Mac, PC, and Linux, as long as the script calling
%   loadNeuroDAP is located inside the NeuroDAP folder or one of its
%   subfolders.
%
%   Example:
%       clear; close all;
%       loadNeuroDAP;
%       [~,~,~,~,~,bluePurpleRed] = loadColors;
%
%   Output:
%       neuroDAPDir - Full path to the detected main NeuroDAP folder.

    %% Find candidate folders for the NeuroDAP root
    %
    % dbstack gives the current call stack.
    % stack(1) is this function: loadNeuroDAP
    % stack(2), if it exists, is the script/function that called loadNeuroDAP.
    stack = dbstack('-completenames');

    candidateStartDirs = {};

    if numel(stack) >= 2
        % Called from a script or function
        callerPath = stack(2).file;
        candidateStartDirs{end+1} = fileparts(callerPath);
    end

    % MATLAB can run editor selections/sections from temporary Editor_*
    % scripts. In that case the caller is outside the repository, so fall
    % back to the folder containing this loader and then the Current Folder.
    candidateStartDirs{end+1} = fileparts(mfilename('fullpath'));
    candidateStartDirs{end+1} = pwd;
    candidateStartDirs = unique(candidateStartDirs, 'stable');

    %% Walk upward from each candidate until we find the NeuroDAP folder
    neuroDAPDir = '';

    for i = 1:numel(candidateStartDirs)
        startDir = candidateStartDirs{i};
        searchDir = startDir;

        while ~isfolder(fullfile(searchDir, 'Methods'))
            parentDir = fileparts(searchDir);

            % If parentDir is the same as searchDir, we reached the filesystem
            % root and did not find a Methods folder.
            if strcmp(parentDir, searchDir)
                searchDir = '';
                break
            end

            searchDir = parentDir;
        end

        if ~isempty(searchDir)
            neuroDAPDir = searchDir;
            break
        end
    end

    if isempty(neuroDAPDir)
        error(['Could not find the NeuroDAP root folder. ', ...
               'Starting folders were: %s'], strjoin(candidateStartDirs, ', '));
    end

    %% Add NeuroDAP/Methods to MATLAB path
    %
    % fullfile automatically handles Mac/Linux "/" versus Windows "\".
    methodsDir = fullfile(neuroDAPDir, 'Methods');

    if ~isfolder(methodsDir)
        error('Found NeuroDAP root, but Methods folder does not exist: %s', methodsDir);
    end

    addpath(genpath(methodsDir));

    %% Optional display message
    fprintf('NeuroDAP loaded from: %s\n', neuroDAPDir);
    fprintf('Added Methods folder to path: %s\n', methodsDir);

end
