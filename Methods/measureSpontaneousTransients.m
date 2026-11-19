function [animalSummary,sessionSummary,windowSummary] = measureSpontaneousTransients(animals,cfg)
% Apply the paper's peak inclusion rule to saved NeuroDAP dF/F snippets.
% https://doi.org/10.64898/2026.07.24.740623 (Fiber photometry methods)
% Raw preprocessing in the paper: 50 Hz lowpass, downsample to 50 Hz,
% sliding 5 ms average, recording-mean normalization. These inputs already
% contain processed dF/F, so we retain their preprocessing. Re-normalizing
% snippets would change the meaning of the configured percent threshold.
%
% Default: Water-aligned Random trials, strictly >3 s after reward and
% <2 s before the next cue. This conservative subset also avoids cues and
% duplicate time across adjacent trials. ITI is stored in seconds. For
% non-reward anchors, caller must explicitly confirm event-free windows.
% Source labels select snippets; measured command activity assigns ON/OFF.
% This measures actuator output, not whether a zero-output controller is armed.
records = cell(0,13);
for selection = 1:numel(cfg.spontaneousSelections)
    sel = cfg.spontaneousSelections(selection);
    assert(sel.window(1)<sel.window(2),'Invalid spontaneous window.');
    rows = contains({animals.task},sel.task,'IgnoreCase',true) & ...
        strcmpi({animals.name},sel.signal) & strcmpi({animals.event},sel.event);
    if ~any(strcmpi(string(cfg.animals),'All'))
        rows = rows & ismember(lower(string({animals.animal})),lower(string(cfg.animals)));
    end
    assert(sel.rewardAligned || sel.confirmedEventFree, ...
        'Confirm custom windows are free of rewards/cues before analyzing.');
    if sel.rewardAligned
        assert(startsWith(sel.event,'Water','IgnoreCase',true), ...
            'Automatic reward exclusion requires a Water-aligned event.');
        assert(sel.window(1)>=cfg.rewardExclusionSeconds, ...
            'Start the window at or after the reward exclusion interval.');
        assert(isfield(sel,'cueToRewardSeconds') && isscalar(sel.cueToRewardSeconds) && ...
            isfinite(sel.cueToRewardSeconds) && sel.cueToRewardSeconds>=0, ...
            'Specify the cue-to-reward delay in seconds for this selection.');
    end
    for r = reshape(find(rows),1,[])
        row = animals(r); t = double(row.timestamp(:)');
        assert(numel(t)>2 && size(row.data,2)==numel(t),'Invalid trace dimensions.');
        dt = median(diff(t));
        assert(dt>0 && max(abs(diff(t)-dt))<dt*0.01,'Irregular timestamps.');
        assert(t(1)<=sel.window(1)+dt/2 && t(end)>=sel.window(2)-dt/2, ...
            'Requested window unavailable for row %d.',r);
        tt = row.trialInfo.trialTable;
        assert(height(tt)==size(row.data,1),'Trial metadata length mismatch.');
        assert(ismember('performing',tt.Properties.VariableNames),'Missing performing flag.');
        if sel.rewardAligned
            assert(all(ismember({'ITI','isReward'},tt.Properties.VariableNames)), ...
                'Water windows require ITI (seconds) and isReward metadata.');
        end
        [clampActive,blueRow,redRow] = getTrialClampActivity(animals,r, ...
            cfg.clampActiveThresholdPercent);
        % Trial-number resets identify concatenated sessions. Keep each
        % row's session labels distinct; selections should use one anchor
        % event per condition, avoiding duplicate representations of trials.
        trialNumbers = row.trialInfo.trialNumber(:);
        sessionIdx = cumsum([true;diff(trialNumbers)<0]);
        sessionNames = string(row.options.sessionList);
        if numel(sessionNames)==1
            sessionIdx = ones(size(trialNumbers));
        elseif numel(sessionNames)~=max(sessionIdx)
            % Some legacy ON/OFF rows map repeated clamp epochs to the
            % nearest behavioral trial, producing false trial-number resets.
            % Preserve the row as an explicitly labeled pooled stratum
            % rather than inventing session boundaries.
            warning('ClampQuality:PooledSessions', ...
                'Row %d: ambiguous legacy session boundaries; pooling this row.',r);
            sessionNames = "Pooled row " + r + ": " + strjoin(sessionNames,', ');
            sessionIdx = ones(size(trialNumbers));
        end
        for j = reshape(find(tt.performing==1),1,[])
            last = sel.window(2);
            if sel.rewardAligned
                if tt.isReward(j)~=1 || ~isfinite(tt.ITI(j)); continue; end
                % End before the next cue, even when it is not a reward.
                last = min(last,tt.ITI(j)-sel.cueToRewardSeconds-cfg.nextCueGuardSeconds);
            end
            idx = t>sel.window(1)+dt*1e-6 & t<last-dt*1e-6;
            if sum(idx)<3; continue; end
            key = sessionNames(sessionIdx(j));
            for state = 0:1
                % Never bridge missing samples or a change in command state.
                valid = idx & clampActive(j,:)==state & isfinite(row.data(j,:));
                edges = diff([false,valid,false]);
                starts = find(edges==1); stops = find(edges==-1)-1;
                for segment = 1:numel(starts)
                    samples = starts(segment):stops(segment);
                    if numel(samples)<3; continue; end
                    x = cfg.scale*double(row.data(j,samples));
                    % 5 ms moving mean is one sample for saved traces <=100 Hz.
                    x = movmean(x,max(1,round(0.005/dt)));
                    amps = detectPaperPositivePeaks(x,cfg.minimumRisePercent);
                    records(end+1,:) = {string(row.animal),state+1,key,r,j, ...
                        numel(amps),numel(x)*dt,sum(amps),mean(amps), ...
                        t(samples(1)),t(samples(end))+dt,blueRow,redRow}; %#ok<AGROW>
                end
            end
        end
    end
end
windowSummary = cell2table(records,'VariableNames',{'Animal','Condition','Session', ...
    'AnimalRow','TrialRow','Count','Seconds','AmplitudeSum','MeanAmplitude', ...
    'StartTime','EndTime','BlueClampRow','RedClampRow'});
windowSummary.Animal = string(windowSummary.Animal);
windowSummary.Session = string(windowSummary.Session);
% One metric per session, as in the paper; ambiguous legacy rows are
% explicitly labeled pooled strata. Animal is the plotting unit.
sessions = cell(0,7);
if ~isempty(windowSummary)
    [g,an,co,se] = findgroups(windowSummary.Animal,windowSummary.Condition,windowSummary.Session);
    for k = 1:numel(an)
        w = windowSummary(g==k,:); count = sum(w.Count); seconds = sum(w.Seconds);
        amp = NaN; if count>0; amp=sum(w.AmplitudeSum)/count; end
        sessions(end+1,:) = {an(k),co(k),se(k),amp,count/seconds,count,seconds}; %#ok<AGROW>
    end
end
sessionSummary = cell2table(sessions,'VariableNames', ...
    {'Animal','Condition','Session','Peak','Hz','Count','Seconds'});
sessionSummary.Animal = string(sessionSummary.Animal);
names = unique(sessionSummary.Animal); peak=nan(numel(names),2); hz=peak;
count=zeros(numel(names),2); seconds=count;
for a=1:numel(names)
    for c=1:2
        rows=sessionSummary.Animal==names(a) & sessionSummary.Condition==c;
        if ~any(rows); continue; end
        peak(a,c)=mean(sessionSummary.Peak(rows),'omitnan');
        hz(a,c)=mean(sessionSummary.Hz(rows));
        count(a,c)=sum(sessionSummary.Count(rows)); seconds(a,c)=sum(sessionSummary.Seconds(rows));
    end
end
animalSummary=table(names,peak(:,1),peak(:,2),hz(:,1),hz(:,2), ...
    count(:,1),count(:,2),seconds(:,1),seconds(:,2), ...
    repmat(cfg.minimumRisePercent,numel(names),1),'VariableNames', ...
    {'Animal','OffPeak','OnPeak','OffHz','OnHz','OffCount','OnCount', ...
    'OffSeconds','OnSeconds','ThresholdPercent'});
end
