function analyzeSessions_Shijia(sessionpath,options)

arguments
    sessionpath string
    options.analyzeTraces logical = true
    options.redo logical = true % Recalculate trial table and all preprocessing
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotLicks logical = true % Plot lick raster summary plot

    options.lick_binSize double = 0.1
end

%% Notes
% Modified from Shun_analyzeBehavior_optoPair
% Shun Li, 11/20/2023
% 02/14/2023: tidied up code, renamed to analyzeBehavior_optoPair
% 2023/07/28: packaged trial table into a function
% 2023/09/02: added camera plotting
% 2023/09/05: changed baselineIdx to selecting baseline licks 
% 2023/10/23: changed how to plot photometry signal, assume everything
% recorded in labjack

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
date = dirsplit{1}; animal = dirsplit{2};
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


%% Save photometry/lick/eye PSTHs

waterIdx = find(labjack.solenoid);
cueIdx = find(labjack.cue);
lickIdx = find(labjack.lick);

if options.analyzeTraces
    analysisEvents = {waterIdx,cueIdx,lickIdx};
    analysisLabels = {'Water','Cue','Lick'};

    analysis = analyzeTraces(timeSeries,labjack.lick,analysisEvents,analysisLabels,params);
end


%% Test plotting traces

timeRange = [-0.5,3];

% Find the number of photometry channels
photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
% photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
nSignals = length(photometryIdx);
disp(['Finished: found ', num2str(nSignals),' photometry signals']);

if options.plotPhotometry

    for photometry = 1:nSignals
        % Load signal of interest
        path = photometryIdx(photometry);
        signal = timeSeries(path).data;
        finalFs = timeSeries(path).finalFs;
        system = timeSeries(path).system;

        initializeFig(.5,.5); tiledlayout('flow');
        [~,~] = plotTraces(waterIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,...
                                signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(cueIdx,timeRange,signal,bluePurpleRed(2,:),params,...
                        signalFs=finalFs,...
                        signalSystem=system,eventSystem=params.session.baselineSystem);
        [~,~] = plotTraces(lickIdx,timeRange,signal,bluePurpleRed(3,:),params,...
                signalFs=finalFs,...
                signalSystem=system,eventSystem=params.session.baselineSystem);

        plotLicks(waterIdx,timeRange,options.lick_binSize,bluePurpleRed(1,:),[],labjack.lick,params);
    end

end
% signalSystem is LJ (or PMT), eventSystem is LJ (or NIDQ)