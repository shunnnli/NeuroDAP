function [glm_ds] = generateGLMTable_forEmilyLJ(sessionName,timeSeries,leftTone,rightTone,rightLick,rightSolenoid,params,trials)
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

%2024/7/11
%Adapted by EF to work with new rig using labjack to collect behavior data


%% Generate table for GLM for each file
% Initialize table value
varTypes = {'string','double','double','double','double','double','double','string','double','double'};
varNames = {'SessionName','TrialNumber','Timestamp','go','nogo','lick','reward','outcome','photometryLhemi','photometryRhemi'};

% Define target frequency and downsample factor
targetFs = params.sync.photometryFs(1);
BehSamplesPerBin = params.sync.behaviorFs / targetFs;
LJSamplesPerBin = params.sync.labjackFs / targetFs;

% Align to first common sync pulse
timePhotometry = params.sync.timePhotometry(params.sync.commonStartPhotometry:end);
commonStartTimeSeries = ceil(params.sync.commonStartPhotometry/LJSamplesPerBin);

% Find common samples
total_ds_Behsamples = floor(length(timePhotometry)/BehSamplesPerBin);
total_ds_TSsamples = length(timeSeries(1).data)-commonStartTimeSeries;
total_ds_samples = min(total_ds_Behsamples, total_ds_TSsamples);
glm_ds = table('Size',[total_ds_samples length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Downsample behavior data & match nSamples
BehWindow = params.sync.commonStartPhotometry : params.sync.commonStartPhotometry+total_ds_samples*BehSamplesPerBin-1;
go_ds = squeeze(sum(...
    reshape(leftTone(BehWindow), 1, BehSamplesPerBin, total_ds_samples), ...
    2))>0;
nogo_ds = squeeze(sum(...
    reshape(rightTone(BehWindow), 1, BehSamplesPerBin, total_ds_samples), ...
    2))>0;
lick_ds = squeeze(sum(...
    reshape(rightLick(BehWindow), 1, BehSamplesPerBin, total_ds_samples), ...
    2))>0;
reward_ds = squeeze(sum(...
    reshape(rightSolenoid(BehWindow), 1, BehSamplesPerBin, total_ds_samples), ...
    2))>0;

%% Add data to table
glm_ds.SessionName(:) = sessionName;
glm_ds.Timestamp(:) = 0:total_ds_samples-1;
glm_ds.go(:) = go_ds;
glm_ds.nogo(:) = nogo_ds;
glm_ds.lick(:) = lick_ds;
glm_ds.reward(:) = reward_ds;

timeSeriesWindow = commonStartTimeSeries:commonStartTimeSeries+total_ds_samples-1;

numTimeSeries = size(timeSeries,2);

if numTimeSeries == 2
    glm_ds.photometryLhemi(:) = timeSeries(1).data(timeSeriesWindow);
    glm_ds.photometryRhemi(:) = timeSeries(2).data(timeSeriesWindow);
else
    %for 20231115 data or data when there is only one hemisphere of data
    % - may need to specify below which hemisphere is being collected
    glm_ds.photometryLhemi(:) = timeSeries.data(timeSeriesWindow); 
end

goCue_ds = find(go_ds); nogoCue_ds = find(nogo_ds);
allTrials_ds = sortrows([[goCue_ds;nogoCue_ds], [ones(length(goCue_ds),1); zeros(length(nogoCue_ds),1)]],1);

%% 1/4/24 - EF edit 

% Clean trials table up (remove last row with no data)
for i = 1:height(trials)
    if trials.TrialNumber(i) == 0
        trials(i,:) = []; 
    end 
end 

% For edge case when commonStartTime is after the first trial or the time 
% of the last cue is after the end of BehWindow (meaning that the number of original trials 
% (in "trials" table)is greater than the the number of downsampled trials (allTrials_ds))

if height(trials) > height(allTrials_ds)
    newFirstTrial = find(trials.CueTime > params.sync.commonStartPhotometry,1);
    newLastTrial = find(trials.CueTime < BehWindow(end));
    newLastTrial = newLastTrial(end); 
    trials = trials(newFirstTrial:newLastTrial,:);
end

%% Fill in trial outcomes for GLM table

for i = 1:height(trials)
    cur_cue = allTrials_ds(i);

    if i == 1; glm_ds.outcome(1:cur_cue) = 'NA'; end
    if i == height(trials); next_cue = total_ds_samples;
    else; next_cue = allTrials_ds(i+1);
    end

    glm_ds.TrialNumber(cur_cue:next_cue) = i;
    glm_ds.outcome(cur_cue:next_cue) = trials{i,'outcome'};
end

% plotTraces(find(glm_ds.go),[-1,3],glm_ds.photometryLhemi,signalFs=50,sameSystem=true,...
%         color=[0.532, 0.2323,0.44],plotIndividual=true);

end

%% For generating GLM table
% Final format: https://github.com/bernardosabatinilab/sabatinilab-glm/blob/refactor/notebooks/glm_fit.ipynb

