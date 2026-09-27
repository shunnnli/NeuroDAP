function [data,info] = rescaleClampTraces(data,options)

% Rescale stored blueClamp/redClamp traces so that the command level the
% clamp actually reaches maps onto 100%.
%
% Sessions loaded before loadSessions auto-detected the clamp full scale
% (autoClampMax) were normalized by a nominal ADC range that can differ from
% the hardware. A red command topping out at 500 counts read against
% [25,600] peaks at 82.6% instead of 100%. Since the clamp reaches the same
% maximum at some point in every session, that peak IS the full scale, so
% the stored percentages can be corrected offline without reloading:
%
%   pct_new = pct_old * 100 / detectedMax
%
% Takes either a summary or an animals struct: rows are grouped by whichever
% of animal/date/session/task/name they carry, and the max is pooled over
% all rows (ie all events) of the same session and channel. Traces and the
% stage statistics derived from them are scaled together.
%
% Rescaling is exact only where the stored trace was not clipped. Values
% above the nominal max were clipped to 100% when the trace was created and
% cannot be recovered - those groups are reported as 'clipped'.
%
% Example:
%   [summary,info] = rescaleClampTraces(summary,dateBefore=20260927);
%   animals = rescaleClampTraces(animals);

arguments
    data struct % summary or animals struct
    options.clampSignals cell = {'redClamp','blueClamp'}
    options.dateBefore double = [] % only rescale sessions before this date (yyyymmdd)
    options.targetPct double = 100 % detected max is mapped onto this percentage
    options.minSamples double = 5  % samples required at the detected max
    options.tolerance double = 0.5 % (%) how close to the max is the same level
    options.minPct double = 10     % skip groups peaking below this (clamp likely off)
    options.verbose logical = true
end

info = table('Size',[0 6],'VariableTypes',{'string','string','double','double','double','string'},...
             'VariableNames',{'group','name','nRows','detectedMax','scaleFactor','status'});
if isempty(data); return; end

% Rows holding a clamp command
signalNames = string({data.name});
clampRows = find(ismember(lower(signalNames),lower(string(options.clampSignals))));
if isempty(clampRows)
    if options.verbose; disp('Finished: no redClamp/blueClamp rows found, nothing to rescale'); end
    return
end

% Only rescale sessions recorded before dateBefore (needs a date field)
if ~isempty(options.dateBefore)
    if isfield(data,'date')
        rowDates = str2double(string({data(clampRows).date}));
        clampRows = clampRows(rowDates < options.dateBefore);
        if isempty(clampRows)
            if options.verbose; disp('Finished: no clamp rows before dateBefore, nothing to rescale'); end
            return
        end
    else
        warning('rescaleClampTraces:NoDateField',...
            ['No date field (animals struct?), dateBefore ignored and all ',...
             num2str(length(clampRows)),' clamp rows are rescaled.']);
    end
end

% Group rows by session & channel, using whichever fields this struct has
keyFields = {'animal','date','session','task','name'};
keyFields = keyFields(isfield(data,keyFields));
groupLabels = strings(length(clampRows),1);
for i = 1:length(clampRows)
    groupLabels(i) = strjoin(string(cellfun(@(f) string(data(clampRows(i)).(f)),...
                                            keyFields,UniformOutput=false)),'_');
end
[uniqueGroups,~,groupIdx] = unique(groupLabels);

for g = 1:length(uniqueGroups)
    groupRows = clampRows(groupIdx == g);
    groupName = string(data(groupRows(1)).name);

    % Pool every trial of this session & channel to find the full scale
    pooled = [];
    for r = reshape(groupRows,1,[])
        if ~isnumeric(data(r).data); continue; end
        pooled = [pooled; data(r).data(:)]; %#ok<AGROW>
    end
    detectedMax = detectPlateauMax(pooled,minSamples=options.minSamples,...
                                   tolerance=options.tolerance);

    % Decide what to do with this group
    if isnan(detectedMax) || detectedMax < options.minPct
        status = "skipped: peak below minPct";
        scaleFactor = 1;
    elseif detectedMax >= options.targetPct
        % Already at (or clipped to) full scale, nothing to recover
        status = "unchanged: already at targetPct";
        scaleFactor = 1;
    else
        scaleFactor = options.targetPct / detectedMax;
        status = "rescaled";
    end

    % Apply to traces and to the stage statistics derived from them
    if scaleFactor ~= 1
        for r = reshape(groupRows,1,[])
            if isnumeric(data(r).data); data(r).data = data(r).data * scaleFactor; end

            rowFields = fieldnames(data(r));
            stageFields = rowFields(contains(rowFields,'stage',IgnoreCase=true));
            for f = 1:length(stageFields)
                cur_stage = data(r).(stageFields{f});
                if isstruct(cur_stage) && isfield(cur_stage,'data') && isnumeric(cur_stage.data)
                    cur_stage.data = cur_stage.data * scaleFactor;
                    data(r).(stageFields{f}) = cur_stage;
                end
            end
        end
    end

    info = [info; {uniqueGroups(g),groupName,length(groupRows),detectedMax,scaleFactor,status}]; %#ok<AGROW>
    if options.verbose
        disp(['Finished: ',char(uniqueGroups(g)),' -> ',char(status),...
              ' (max ',num2str(detectedMax,'%.2f'),'%, x',num2str(scaleFactor,'%.4f'),')']);
    end
end

if options.verbose
    nRescaled = sum(info.status == "rescaled");
    disp(['Finished: rescaled ',num2str(nRescaled),'/',num2str(height(info)),...
          ' clamp session-channels to ',num2str(options.targetPct),'% full scale']);
end

end
