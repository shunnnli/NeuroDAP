function protocol = getCellProtocol(headerString,options)

% Shun Li, 12/21/2023
% Modified from phUtil_rcCheck_Rin by Bernardo Sabatini
% Get QC (Rs, Rm, Cm) for each sweep

arguments
    headerString

    options.outputFs double = 10000
    options.rcCheckRecoveryWindow double = 100 % in ms
end

% Get cycle
protocol.cycle = getHeaderValue(headerString,'state.cycle.cycleName');

% Get experiment protocol
% UNFINISHED: ADD BLUE AND RED SEPARATION
protocol.amplitude = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='amplitude');
protocol.duration = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='duration');
protocol.offset = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='offset');
protocol.numPulses = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='numPulses');
protocol.isi = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='isi');
protocol.pulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='pulseWidth');
protocol.delay = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='delay');
protocol.ramp = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='ramp');
protocol.patternRepeats = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='patternRepeats');
protocol.patternISI = getHeaderValue(headerString,'state.phys.internal.pulseString_ao1',pulseVar='patternISI');

% Calculate stim onset time
protocol.stimOnset = protocol.delay + ((1:protocol.numPulses)-1)*protocol.isi;
protocol.stimOnset = protocol.stimOnset * (options.outputFs/1000); % in samples

% Extract DMD related params
if contains(protocol.cycle,'randomSearch')
    protocol.depth = getHeaderValue(headerString,'state.zDMD.searchDepth',convert=true);
    protocol.repetition = getHeaderValue(headerString,'state.zDMD.searchRepetition',convert=true);
end

end