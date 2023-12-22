function [glm_ds] = generateGLMTable_forEmily(sessionName,timeSeries,leftTone,rightTone,rightLick,rightSolenoid,params,trials)
arguments
    sessionName string
    timeSeries struct
    leftTone double
    rightTone double
    rightLick logical
    rightSolenoid double
    params struct
    trials table
end

% Shun Li, 2023/07/25
% Code for emily to generate tables for Josh's GLM analysis for each
% selected sections

% 2023/7/25
% only finished the NI photometry part, LJ photometry part have sample rate
% mismatch issues

%2023/12/07
%updated for compatibility with most recent photometry params + for LJ
%photometry

%2023/12/10
%Adapted as function for NeuroDAP


%% Generate table for GLM for each file
% Initialize table value
varTypes = {'string','double','double','double','double','double','double','string','double','double'};
varNames = {'SessionName','TrialNumber','Timestamp','go','nogo','lick','reward','outcome','photometryLhemi','photometryRhemi'};

% Define target frequency and downsample factor
targetFs = params.sync.photometryFs(1);
NISamplesPerBin = params.sync.behaviorFs / targetFs;
LJSamplesPerBin = params.sync.labjackFs / targetFs;

% Align to first common sync pulse
timeNI = params.sync.timeNI(params.sync.commonStartNI:end);
commonStartTimeSeries = ceil(params.sync.commonStartPhotometry/LJSamplesPerBin);

% Find common samples
total_ds_NIsamples = floor(length(timeNI)/NISamplesPerBin);
total_ds_TSsamples = length(timeSeries(1).data)-commonStartTimeSeries;
total_ds_samples = min(total_ds_NIsamples, total_ds_TSsamples);
glm_ds = table('Size',[total_ds_samples length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Downsample behavior data & match nSamples
go_ds = squeeze(sum(...
    reshape(leftTone(1:total_ds_samples*NISamplesPerBin), 1, NISamplesPerBin, total_ds_samples), ...
    2))>0;
nogo_ds = squeeze(sum(...
    reshape(rightTone(1:total_ds_samples*NISamplesPerBin), 1, NISamplesPerBin, total_ds_samples), ...
    2))>0;
lick_ds = squeeze(sum(...
    reshape(rightLick(1:total_ds_samples*NISamplesPerBin), 1, NISamplesPerBin, total_ds_samples), ...
    2))>0;
reward_ds = squeeze(sum(...
    reshape(rightSolenoid(1:total_ds_samples*NISamplesPerBin), 1, NISamplesPerBin, total_ds_samples), ...
    2))>0;

% Add data to table
glm_ds.SessionName(:) = sessionName;
glm_ds.Timestamp(:) = 0:total_ds_samples-1;
glm_ds.go(:) = go_ds;
glm_ds.nogo(:) = nogo_ds;
glm_ds.lick(:) = lick_ds;
glm_ds.reward(:) = reward_ds;

timeSeriesWindow = commonStartTimeSeries:commonStartTimeSeries+total_ds_samples-1;
glm_ds.photometryLhemi(:) = timeSeries(1).data(timeSeriesWindow);
glm_ds.photometryRhemi(:) = timeSeries(2).data(timeSeriesWindow);

goCue_ds = find(go_ds); nogoCue_ds = find(nogo_ds);
allTrials_ds = sortrows([[goCue_ds;nogoCue_ds], [ones(length(goCue_ds),1); zeros(length(nogoCue_ds),1)]],1);

% Filling in TrialNumber, choice, outcome
for i = 1:height(trials)
    cur_cue = allTrials_ds(i);

    if i == 1; glm_ds.outcome(1:cur_cue) = 'NA'; end
    if i == height(trials); next_cue = total_ds_samples;
    else; next_cue = allTrials_ds(i+1);
    end

    glm_ds.TrialNumber(cur_cue:next_cue) = i;
    glm_ds.outcome(cur_cue:next_cue) = trials{i,'outcome'};
end
    
end

%% For generating GLM table
% Final format: https://github.com/bernardosabatinilab/sabatinilab-glm/blob/refactor/notebooks/glm_fit.ipynb

