function [events, onsetByCombo] = extractPWMStim(clampTrace, options)
% Segment stims by rising edges (> interStimMin_s gap) then compute duty from RAW analog.
% 'duty' in the output table is in **PWM COUNTS** (0..2^pwmResolutionBits-1).

arguments
    clampTrace double
    options.fs (1,1) double         = 10000
    options.label (1,1) string      = "blue"
    options.duration double         = []      % e.g., [0.1 0.5 1]
    options.interStimMin_s double   = 1.0
    options.lowPct double           = 0.1
    options.highPct double          = 99.999999999999999999999
    options.pwmResolutionBits double = 8      % 8-bit → 0..255
    options.pwm double              = []      % e.g., [20 30 40] to snap to known settings
end

Fs    = options.fs;
label = options.label;
x     = clampTrace(:);

% ---- Rising edges from the analog trace (make column + use vertical cat) ----
riseVec = (x > (max(x)/2));
temp    = [false; diff(riseVec)];    % NOTE: semicolon for vertical concat
riseVec = (temp == 1);

% -------- 1) Segment stims from rising edges (gap > interStimMin_s) --------
idx = find(riseVec == 1);
if isempty(idx)
    events = table(string.empty(0,1), [], [], [], [], ...
        'VariableNames', {'label','onset','duration_s','duty','power_pct'});
    onsetByCombo = table(string.empty(0,1), cell(0,1), ...
        'VariableNames', {'combo','onsets'});
    return
end
gapSamp = round(options.interStimMin_s * Fs);
grp     = cumsum([true; diff(idx) > gapSamp]);
cells   = accumarray(grp, idx, [], @(v){v});
nG      = numel(cells);

% -------- 2) For each stim: onset/duration from edges, duty from RAW analog --------
V0  = prctile(x, options.lowPct);
V1  = prctile(x, options.highPct);
if V1 <= V0, V1 = max(x); V0 = min(x); end

onset    = zeros(nG,1);
duration_s = zeros(nG,1);
duty_frac  = zeros(nG,1);  % 0..1

for i = 1:nG
    v = cells{i};
    onset_idx = v(1);

    % include the last carrier cycle
    if numel(v) > 1
        perSamp = round(median(diff(v)));
    else
        perSamp = 0;
    end
    end_idx = min(length(x), v(end) + max(0,perSamp));

    seg = onset_idx:end_idx;
    mu  = mean(x(seg));

    onset(i)    = onset_idx;
    duration_s(i) = numel(seg) / Fs;
    duty_frac(i)  = max(0, min(1, (mu - V0) / (V1 - V0)));
end

% ---- 3) Convert duty to COUNTS, optionally snap to known values ----
maxCount    = 2^options.pwmResolutionBits - 1;    % 255 for 8-bit
duty_counts = round(duty_frac * maxCount);        % e.g., 20/30/40
if ~isempty(options.pwm)
    [~,ix] = min(abs(duty_counts - options.pwm(:)'), [], 2);
    duty_counts = options.pwm(ix);
end
power_pct = 100 * duty_frac;

% Optional: snap duration to canonical set
if ~isempty(options.duration)
    [~,ix]   = min(abs(duration_s - options.duration(:)'), [], 2);
    duration_s = options.duration(ix);
end

% ---- 4) Build outputs (ensure columns, equal lengths) ----
onset     = onset(:);
duration_s  = duration_s(:);
duty_counts = duty_counts(:);
power_pct   = power_pct(:);

n = numel(onset);
lbl = repmat(label, n, 1);

events = table(lbl, onset, duration_s, duty_counts, power_pct, ...
    'VariableNames', {'label','onset','duration_s','duty','power_pct'});
events = sortrows(events,'onset');

% -------- 5) Group onsets by (label, power bin, duration) --------
pwrKey = round(power_pct/5)*5;        % 0,5,10,... (%)
durKey = round(duration_s, 3);
combo  = compose("%s_%dpct_%0.3fs", lbl, pwrKey, durKey);

[G, keys] = findgroups(combo);
onsets = splitapply(@(x){x}, events.onset(:), G);
onsetByCombo = table(keys, onsets, 'VariableNames', {'combo','onsets_s'});
end
