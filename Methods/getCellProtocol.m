function protocol = getCellProtocol(headerString,options)

% Shun Li, 12/21/2023
% Modified from phUtil_rcCheck_Rin by Bernardo Sabatini
% Get QC (Rs, Rm, Cm) for each sweep

arguments
    headerString

    options.outputFs double = 10000
    options.rcCheckRecoveryWindow double = 100 % in ms

    options.rig string
end

% Get cycle
protocol.cycle = getHeaderValue(headerString,'state.cycle.cycleName');

% Decide color
protocol.lastLinesUsed = eval(getHeaderValue(headerString,'state.phys.internal.lastLinesUsed'));
protocol.lastPulsesUsed = eval(getHeaderValue(headerString,'state.phys.internal.lastPulsesUsed'));
protocol.stimChannel = protocol.lastLinesUsed{end};

% Get experiment protocol
pulseString = ['state.phys.internal.pulseString_',protocol.stimChannel];
protocol.amplitude = getHeaderValue(headerString,pulseString,pulseVar='amplitude');
protocol.duration = getHeaderValue(headerString,pulseString,pulseVar='duration');
protocol.offset = getHeaderValue(headerString,pulseString,pulseVar='offset');
protocol.numPulses = getHeaderValue(headerString,pulseString,pulseVar='numPulses');
protocol.isi = getHeaderValue(headerString,pulseString,pulseVar='isi');
protocol.pulseWidth = getHeaderValue(headerString,pulseString,pulseVar='pulseWidth');
protocol.delay = getHeaderValue(headerString,pulseString,pulseVar='delay');
protocol.ramp = getHeaderValue(headerString,pulseString,pulseVar='ramp');
protocol.patternRepeats = getHeaderValue(headerString,pulseString,pulseVar='patternRepeats');
protocol.patternISI = getHeaderValue(headerString,pulseString,pulseVar='patternISI');

% Calculate stim onset time
protocol.stimOnset = protocol.delay + ((1:protocol.numPulses)-1)*protocol.isi;
protocol.stimOnset = protocol.stimOnset * (options.outputFs/1000); % in samples

% Get RC check params
if isempty(getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck'))
    if strcmpi(options.rig,'Wengang')
        rcString = 'amplitude=-5;duration=3000;offset=0;numPulses=1;isi=0;pulseWidth=100;delay=2800;ramp=0;patternRepeats=0;patternISI=0;';
        rcCheckOnset = getHeaderValue(rcString,'',pulseVar='delay') * (options.outputFs/1000);
        rcCheckPulseWidth = getHeaderValue(rcString,'',pulseVar='pulseWidth') * (options.outputFs/1000);
        rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
        rcCheckAmplitude = getHeaderValue(rcString,'',pulseVar='amplitude');
        protocol.rcCheckOnset = rcCheckOnset;
        protocol.rcCheckPulseWidth = rcCheckPulseWidth;
        protocol.rcCheckEnd = rcCheckEnd;
        protocol.rcCheckAmplitude = rcCheckAmplitude;
    elseif strcmpi(options.rig,'Paolo')
        rcString = 'amplitude=-10;duration=1000;offset=0;numPulses=1;isi=0;pulseWidth=20;delay=100;ramp=0;patternRepeats=0;patternISI=0;';
        rcCheckOnset = getHeaderValue(rcString,'',pulseVar='delay') * (options.outputFs/1000);
        rcCheckPulseWidth = getHeaderValue(rcString,'',pulseVar='pulseWidth') * (options.outputFs/1000);
        rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
        rcCheckAmplitude = getHeaderValue(rcString,'',pulseVar='amplitude');
        protocol.rcCheckOnset = rcCheckOnset;
        protocol.rcCheckPulseWidth = rcCheckPulseWidth;
        protocol.rcCheckEnd = rcCheckEnd;
        protocol.rcCheckAmplitude = rcCheckAmplitude;
    end
else
    rcCheckOnset = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='delay') * (options.outputFs/1000);
    rcCheckPulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='pulseWidth') * (options.outputFs/1000);
    rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
    rcCheckAmplitude = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='amplitude');
    protocol.rcCheckOnset = rcCheckOnset;
    protocol.rcCheckPulseWidth = rcCheckPulseWidth;
    protocol.rcCheckEnd = rcCheckEnd;
    protocol.rcCheckAmplitude = rcCheckAmplitude;
end

% Extract DMD related params
if contains(protocol.cycle,'randomSearch')
    protocol.depth = getHeaderValue(headerString,'state.zDMD.searchDepth',convert=true);
    protocol.repetition = getHeaderValue(headerString,'state.zDMD.searchRepetition',convert=true);

    % Extract cell location
    if strcmp(protocol.stimChannel,'ao1')
        tVector = eval(getHeaderValue(headerString,'state.zDMD.tVectorBlue'));
        protocol.cellX = 342 + round(tVector(1));
        protocol.cellY = 304 + round(tVector(2));
    elseif strcmp(protocol.stimChannel,'ao2')
        tVector = eval(getHeaderValue(headerString,'state.zDMD.tVectorRed'));
        protocol.cellX = 342 - round(tVector(2));
        protocol.cellY = 304 + round(tVector(1));
    end
end

% Save other values
protocol.headerString = headerString;
protocol.pulseString = pulseString;

end