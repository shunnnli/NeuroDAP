function event_wvf = getLFPWaveforms(raw,timeRange,params,options)

arguments
    raw
    timeRange double
    params struct
    options.event double
    options.Ch = 1:385 % range of channel
end


lfpFs = params.sync.lfpFs;
timeLFP = params.sync.timeLFP;
timeNI = params.sync.timeNI;
timesteps = floor(lfpFs*(timeRange(2)-timeRange(1)));

% wvf = number of channels x number of optostim x timesteps
if isempty(options.event)
    nTrials = 500;
    event_wvf = nan(length(options.Ch),nTrials,timesteps);
else
    nTrials = length(options.event);
    event_wvf = nan(length(options.Ch),nTrials,timesteps);
end

bar = waitbar(0,'getLFPWaveforms...','Name','Extracting LFP waveforms',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(bar,'canceling',0);

try

    if ~isempty(options.event)
        % Get opto-triggered waveform
        for j = 1:nTrials
            % Check for clicked Cancel button
            if getappdata(bar,'canceling'); break; end
            % Update waitbar and message
            progress = j/nTrials;
            waitbar(progress,bar,['Event ',num2str(j),'/',num2str(nTrials)]);
    
            % Find Imec stim index
            niStimIdx = options.event(j);
            [~, lfpStimIdx] = min(abs(timeLFP-timeNI(niStimIdx)));
    
            lfpFirstIdx = lfpStimIdx + floor(lfpFs*timeRange(1));
            lfpLastIdx = lfpFirstIdx + timesteps - 1;
        
            wf = raw(options.Ch, lfpFirstIdx:lfpLastIdx);
            % Spike waveform for each stim
            event_wvf(:,j,:) = wf;
        end
    else
        for j = 1:nTrials
            % Check for clicked Cancel button
            if getappdata(bar,'canceling'); break; end
            % Update waitbar and message
            progress = j/nTrials;
            waitbar(progress,bar,['Event ',num2str(j),'/',num2str(nTrials)]);
    
            % Find Imec stim index
            lfpStimIdx = randi(length(timeLFP),1);
    
            lfpFirstIdx = lfpStimIdx + floor(lfpFs*timeRange(1));
            lfpLastIdx = lfpFirstIdx + timesteps - 1;
        
            wf = raw(options.Ch, lfpFirstIdx:lfpLastIdx);
            % Spike waveform for each stim
            event_wvf(:,j,:) = wf;
        end
    end

    delete(bar);

catch ME
    delete(bar);
    rethrow(ME) %re-issue the error
end

% Convert to mV using conversion factor
% event_wvf = event_wvf .* 2.3438;

end