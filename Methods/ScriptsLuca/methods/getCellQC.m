function qc = getCellQC(data, options)

% Shun Li, 12/21/2023
% Modified from phUtil_rcCheck_Rin by Bernardo Sabatini
% Get QC (Rs, Rm, Cm) for each sweep

arguments
    data double

    options.calculate logical = true % if false, just extract relevant value from headerString

    options.headerString
    options.onset double
    options.pulseWidth double

    options.baseline double = 100 % in ms before RC pulse
    options.Fs double = 10000
    options.currentClamp logical = false
    options.amplitude double = -5
    options.mode char {mustBeMember(options.mode,['peak', 'last10%', 'last30%', 'last50%'])} = 'last30%'
end

%% Check input

if ~isfield(options,'headerString')
    if ~isfield(options,["onset","pulseWidth"])
        error('Have to provide onset and pulseWith if headerString not provided!');
    end
    if ~options.calculate
        error("Have to provide headerString if not calculate from raw trace!");
    end
else
    options.currentClamp = getHeaderValue(options.headerString,'state.phys.settings.currentClamp0');
    options.Fs = getHeaderValue(options.headerString,'state.phys.settings.inputRate', convert=true); % ms

    options.amplitude = getHeaderValue(options.headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='amplitude');
    options.onset = getHeaderValue(options.headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='delay');
    options.pulseWidth = getHeaderValue(options.headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='pulseWidth');
end

%% Extract from headerString
if ~options.calculate
    qc.Rs = getHeaderValue(options.headerString,'state.phys.cellParams.rs0',convert=true);
    qc.Rm = getHeaderValue(options.headerString,'state.phys.cellParams.rm0',convert=true);
    qc.Cm = getHeaderValue(options.headerString,'state.phys.cellParams.cm0',convert=true);
    return
end

%% Calculate basic params

time_step = 1000/options.Fs;
start_point = options.onset/time_step+1;
duration = options.pulseWidth/time_step;
end_point = start_point + duration;
end_baseline = max(start_point-1, 1);

start_baseline = (start_point-options.baseline/time_step);
start_baseline=max(start_baseline, 1);

%% Extract baseline
baseline = mean(data(start_baseline:end_baseline));

%% Extract RC check
switch options.mode
    case 'peak'
        if options.amplitude < 0 
            vv=min(data(start_point:end_point));
        else
            vv=max(data(start_point:end_point));
        end
    case 'last10%'
        ss1=start_point+floor((1-0.1)*duration);
        vv=mean(data(ss1:end_point));
    case 'last30%'
        ss1=start_point+floor((1-0.3)*duration);
        vv=mean(data(ss1:end_point));
    case 'last50%'
        ss1=start_point+floor((1-0.5)*duration);
        vv=mean(data(ss1:end_point));
end

%% Calculate Rin (Rm)
if options.currentClamp
    qc.Rm = 1000*(vv-baseline)/options.amplitude;                
else
    qc.Rm = 1000*options.amplitude/(vv-baseline);                
end

end
