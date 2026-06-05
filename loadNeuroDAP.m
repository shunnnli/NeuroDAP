function neuroDAPDir = loadNeuroDAP()
% loadNeuroDAP  Add the NeuroDAP Methods folder to the MATLAB path.
%
%   neuroDAPDir = loadNeuroDAP()
%
%   This function automatically finds the main NeuroDAP folder by starting
%   from the folder of the script/function that called loadNeuroDAP, then
%   walking upward until it finds a folder that contains "Methods".
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

    %% Find the folder of the file that called loadNeuroDAP
    %
    % dbstack gives the current call stack.
    % stack(1) is this function: loadNeuroDAP
    % stack(2), if it exists, is the script/function that called loadNeuroDAP.
    stack = dbstack('-completenames');

    if numel(stack) >= 2
        % Called from a script or function
        callerPath = stack(2).file;
        startDir = fileparts(callerPath);
    else
        % Called directly from the Command Window
        startDir = pwd;
    end

    %% Walk upward until we find the NeuroDAP folder
    %
    % We define the NeuroDAP root as the nearest parent folder containing
    % a subfolder named "Methods".
    neuroDAPDir = startDir;

    while ~isfolder(fullfile(neuroDAPDir, 'Methods'))

        parentDir = fileparts(neuroDAPDir);

        % If parentDir is the same as neuroDAPDir, we reached the filesystem
        % root and did not find a Methods folder.
        if strcmp(parentDir, neuroDAPDir)
            error(['Could not find the NeuroDAP root folder. ', ...
                   'Starting folder was: %s'], startDir);
        end

        neuroDAPDir = parentDir;
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