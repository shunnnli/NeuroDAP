function analyzeSessions_FFT(sessionpath,options)

arguments
    sessionpath string
    options.fftWindow double = 60 % in sec
    options.fillmissing logical = true
end

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; 
if ispc; projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end
% Get animal name and session date
dirsplit = strsplit(sessionName,'-');
date = dirsplit{1}; animal = dirsplit{2}; sessionTask = dirsplit{3};
clear dirsplit

disp(strcat('**********',sessionName,'**********'));
load(strcat(sessionpath,filesep,'timeseries_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'behavior_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'));

if ~isfield(params.session,'name'); params.session.name = sessionName; end
if ~isfield(params.session,'date'); params.session.date = date; end
if ~isfield(params.session,'animal'); params.session.animal = animal; end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end

% Create analysis.mat
if ~isempty(dir(fullfile(sessionpath,"analysis_*.mat")))
    load(strcat(sessionpath,filesep,'analysis_',sessionName,'.mat'));
else
    save(strcat(sessionpath,filesep,'analysis_',sessionName),'sessionName','-v7.3');
    disp('Finished: analysis_.mat not found, created a new one');
end
disp(['Finished: Session ',sessionName,' loaded']);

% Load camera signal
eyeAreaIdx = find(cellfun(@(x) strcmpi(x,'eyeArea'), {timeSeries.name}));
if ~isempty(eyeAreaIdx); eyeArea_detrend = timeSeries(eyeAreaIdx).data; end
pupilAreaIdx = find(cellfun(@(x) strcmpi(x,'pupilArea'), {timeSeries.name}));
if ~isempty(pupilAreaIdx); pupilArea_detrend = timeSeries(pupilAreaIdx).data; end

save(strcat(sessionpath,filesep,'sync_',params.session.name),'params','-append');

%% Perform FFT for the whole recording

initializeFig(0.6,0.5);
nexttile; [~,~] = plotFFT(timeSeries(1).data,Fs=timeSeries(1).finalFs); 
title('NAc'); 
nexttile; [~,~] = plotFFT(timeSeries(2).data,Fs=timeSeries(2).finalFs); 
title('LHb'); 
saveas(gcf,strcat(sessionpath,filesep,'FFT_overall.png'));

%% Perform FFT for rolling time window

initializeFig(.6,.5);
for i = 1:2
    nexttile;
    if options.fillmissing; signal = fillmissing(timeSeries(i).data,'next'); end
    finalFs = timeSeries(i).finalFs;
    spectrogram(signal,options.fftWindow*finalFs,options.fftWindow*finalFs/2,[],finalFs,'yaxis');
    box off
    title(timeSeries(i).name);
end
saveas(gcf,strcat(sessionpath,filesep,'FFT_spectrogram.png'));

end