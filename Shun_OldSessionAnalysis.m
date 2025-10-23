% Shun_oldSessionAnalysis.m
% 2025/04/07

% Analysis old experiments before 2024
% Use old sync.mat file to run analyzeTraces function

%% Load data

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'))';

%% Run through all sessions

% Run each session
errorSessionIdx = []; errorMessage = {};
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx errorMessage taskList

    idx = s; % if loop through all
    sessionpath = sessionList{s};

    dirsplit = strsplit(sessionpath,filesep); 
    sessionName = dirsplit{end}; clear dirsplit

    % Name all output files
    syncOutputName = strcat(sessionpath,filesep,'sync_',sessionName);
    timeseriesOutputName = strcat(sessionpath,filesep,'timeseries_',sessionName);
    dataOutputName = strcat(sessionpath,filesep,'data_',sessionName);
    behaviorOutputName = strcat(sessionpath,filesep,'behavior_',sessionName);

    try
        % Load old sync_.mat from old analysis folder
        load(strcat(sessionList{s},filesep,'Old analysis',filesep,'sync_',sessionName,'.mat'));
        save(behaviorOutputName,'sessionName','session','-v7.3');
        save(syncOutputName,'sessionName','session','-v7.3');
        save(timeseriesOutputName,'sessionName','session','-v7.3');
        save(dataOutputName,'sessionName','session','-v7.3');

        % Build timeSeries
        timeSeries(1).name = 'dLight';
        timeSeries(1).data = demodGreenLP;
        timeSeries(1).finalFs = params.photometry.finalFs;
        timeSeries(1).system = 'LJ';
        timeSeries(1).time_offset = timePhotometry(1) - timeNI(1);
        timeSeries(1).demux = true;
        timeSeries(1).demux_freq = params.photometry.freq1;
        timeSeries(1).detrend = true;
        timeSeries(1).detrend_type = 'rolling-z';
        timeSeries(1).detrend_window = params.photometry.detrendWindowTime;
        timeSeries(1).options = params.photometry;
        save(timeseriesOutputName,'timeSeries','-append');

        % Save nidq stuff to data.mat
        save(dataOutputName,...
            'airpuff','airpuff_rounded','firstPulse','leftLick','rightLick',...
            'leftTone','rightTone','leftTone_rouned','rightTone_rounded','allTones',....
            'leftSolenoid','rightSolenoid',...
            'photometry_raw','blueLaser','redLaser','-append');

        % Save sync data
        save(syncOutputName,'params','-append');


        %% Find event timestamp
        waterIdx = find(rightSolenoid_rounded);
        airpuffIdx = find(airpuff_rounded);

        if strcmp(task,'random')
            toneIdx = find(leftTone); stimIdx = firstPulse;

            randomMinSample = 15*params.sync.behaviorFs;
            randomMaxSample = length(params.sync.timeNI) - (15*params.sync.behaviorFs);
            baselineIdx = sort(randi([randomMinSample,randomMaxSample],100,1));

            stageTime = [-2,0;0,2];
            analysisEvents = {waterIdx,waterLickIdx,airpuffIdx,toneIdx,stimIdx,baselineIdx};
            eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                                findTrials(airpuffIdx,trials),findTrials(toneIdx,trials),...
                                findTrials(stimIdx,trials),findTrials(baselineIdx,trials)};
            
            analysisLabels = {'Water','Rewarded licks','Airpuff','Tone','Stim','Baseline','Blue stim'};
            
            eventTrialNum = eventTrialNum(~cellfun('isempty',analysisEvents));
            analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
            analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));
            
            for i = 1:length(analysisEvents)
                disp(['Total ',analysisLabels{i},': ',num2str(length(analysisEvents{i}))]);
            end
        else
            stimTrials = trials{trials.isTone == 0 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            pairTrials = trials{trials.isTone == 1 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            toneTrials = trials{trials.isTone == 1 & trials.isStim == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
            stimIdx = stimTrials(:,2);
            pairIdx = pairTrials(:,2);
            toneIdx = toneTrials(:,2);

            stageTime = [-2,0;0,0.5;0.5,5];
            analysisEvents = {waterIdx,waterLickIdx,toneIdx,stimIdx,pairIdx,airpuffIdx,baselineIdx};
            eventTrialNum = {findTrials(waterIdx,trials),findTrials(waterLickIdx,trials),...
                            toneTrials(:,1),stimTrials(:,1),pairTrials(:,1),...
                            findTrials(airpuffIdx,trials),findTrials(baselineIdx,trials)};
    
            analysisLabels = {'Water','Rewarded licks','Tone only','Stim only','Pair','Airpuff','Baseline','Blue stim'};
        
            eventTrialNum = eventTrialNum(~cellfun('isempty',analysisEvents));
            analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
            analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));
        
            for i = 1:length(analysisEvents)
                disp(['Total ',analysisLabels{i},': ',num2str(length(analysisEvents{i}))]);
            end
        end

        save(behaviorOutputName,'allTrials','trials',...
            'waterIdx','waterLickIdx','airpuffIdx','toneIdx','stimIdx','pairIdx',...
            'baselineIdx',...
            '-append');

        %% Run analyzeTraces
        
        analysis = analyzeTraces(timeSeries,rightLick,analysisEvents,analysisLabels,params,...
                             stageTime=stageTime,...
                             trialNumber=eventTrialNum,trialTable=trials);

    catch ME
        errorSessionIdx = [errorSessionIdx;s];
        msg = getReport(ME); 
        errorMessage{end+1} = msg; disp(msg);
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end

close all;


%% Load data

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'))';

% LP all 