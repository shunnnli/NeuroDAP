function criteria = dmdDefaultCriteria()

% The hotspot criteria shared by loadSlicesDMD, analyzeDMDSearch and the tests.
%
% Every value mirrors the live configuration of pyNeuroDAP's
% Shun_dmd_analysis notebook, so the two pipelines classify the same spots.
% Keeping them in one place is what stops MATLAB and the notebook drifting
% apart silently: a spot near threshold flips on a change to any of these.

criteria = struct( ...
    'baselineMs',           [-30 -1], ... % levelling and noise, one window
    'thresholdFactor',      3, ...        % x baseline SD
    'minSamples',           50, ...       % 5 ms at 10 kHz, not consecutive
    'reanalysisWindowMs',   [0 50], ...
    'minAmplitudePa',       2, ...        % floor on the window's peak
    'detectionWindowMs',    [5 50], ...   % starts at 5 ms to clear the artifact
    'detectionGapMs',       5, ...        % dead time before the first control
    'derivativeLagMs',      2, ...
    'minChargeZ',           0.5, ...      % on the cell-null z scale
    'minSpotsForReference', 5, ...        % below this a MAD is noise
    'minSweeps',            2, ...        % MIN_SWEEPS_PER_HOTSPOT at QC time
    'peakAverageMs',        1, ...
    'maxRs',                30);          % MOhm, on min(Rs, Rs_headerString)

end
