function [active,blueRow,redRow] = getTrialClampActivity(animals,signalRow,thresholdPercent)
% Command activity on the photometry sample grid: 1=ON, 0=OFF, NaN=unknown.
% Match animal/task/event and verify trial/session identity before combining.
% Commands use saved 0-100 percent units. Either channel above threshold is
% active. Zero command does not distinguish a disabled from an idle controller.
assert(isscalar(thresholdPercent) && isfinite(thresholdPercent) && ...
    thresholdPercent>=0,'Invalid command threshold.');
row = animals(signalRow); t = double(row.timestamp(:)');
matches = strcmpi({animals.animal},row.animal) & ...
    strcmpi({animals.task},row.task) & strcmpi({animals.event},row.event);
commands = cell(1,2); names = {'blueClamp','redClamp'}; commandRows = zeros(1,2);
for channel = 1:2
    found = find(matches & strcmpi({animals.name},names{channel}));
    assert(isscalar(found),'Expected one matching %s row for %s/%s/%s.', ...
        names{channel},row.animal,row.task,row.event);
    command = animals(found); commandRows(channel) = found;
    assert(size(command.data,1)==size(row.data,1) && ...
        isequaln(command.trialInfo.trialNumber(:),row.trialInfo.trialNumber(:)) && ...
        isequal(string(command.options.sessionList(:)),string(row.options.sessionList(:))), ...
        'ClampQuality:CommandTrialMismatch', ...
        'Command trial/session identity mismatch for %s.',names{channel});
    a = command.trialInfo.trialTable; b = row.trialInfo.trialTable;
    assert(height(a)==height(b),'Command trial table length mismatch.');
    if ismember('CueTime',a.Properties.VariableNames) && ismember('CueTime',b.Properties.VariableNames)
        assert(isequaln(a.CueTime,b.CueTime),'Command CueTime mismatch.');
    end
    tc = double(command.timestamp(:)');
    assert(numel(tc)>1 && all(isfinite(tc)) && all(diff(tc)>0) && ...
        size(command.data,2)==numel(tc),'Invalid command timestamps.');
    % Nearest sample preserves command states; no extrapolation or bridging
    % of missing command samples. Commands may have a different sample rate.
    commands{channel} = interp1(tc,double(command.data)',t(:),'nearest',NaN)';
end
blue = commands{1}; red = commands{2};
active = nan(size(blue)); valid = isfinite(blue) & isfinite(red);
active(valid) = blue(valid)>thresholdPercent | red(valid)>thresholdPercent;
blueRow = commandRows(1); redRow = commandRows(2);
end
