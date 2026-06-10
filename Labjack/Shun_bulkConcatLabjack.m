%% Concatenate Labjack files
% Faster locally

clear; close all; loadNeuroDAP;
% addpath(addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\NeuroDAP\Methods')));

% Select sessions
sessionList = uipickfiles('FilterSpec','C:\Shun\Recordings');
errorSessionIdx = []; errorMessage = {};

%% Loop through and concatenate

for i = 1:length(sessionList)
    try
        disp(['Concatenating session: ',num2str(i)]);
        concatLabjack(sessionList{i},save=true,record=[1,1,0],plot=false);
    catch ME
        errorSessionIdx = [errorSessionIdx;i];
        msg = getReport(ME); 
        errorMessage{end+1} = msg; disp(msg);
        warning(['Session ', num2str(i), ' have an error, skipped for now!!!!']);
        continue
    end
end
disp('Finished concatenating session!'); 