function criteria = dmdDefaultCriteria()

% The hotspot criteria shared by loadSlicesDMD, analyzeDMDSearch and the tests.
%
% Follow the notebook's later "Update configuration if needed" cell:
% evoked hotspots, one surviving sweep, and a 5-50 ms metric window.

criteria = struct( ...
    'baselineMs',           [-30 -1], ... % levelling and noise, one window
    'thresholdFactor',      3, ...        % x baseline SD
    'minSamples',           50, ...       % 5 ms at 10 kHz, not consecutive
    'reanalysisWindowMs',   [0 50], ...
    'minAmplitudePa',       2, ...        % floor on the window's peak
    'detectionWindowMs',    [5 50], ...   % starts at 5 ms to clear the artifact
    'metricWindowMs',       [5 50], ...   % net peak and AUC heatmap window
    'detectionGapMs',       5, ...        % dead time before the first control
    'derivativeLagMs',      2, ...
    'minChargeZ',           0.5, ...      % on the cell-null z scale
    'minSpotsForReference', 5, ...        % below this a MAD is noise
    'minSweeps',            1, ...        % later MIN_SWEEPS_PER_HOTSPOT setting
    'peakAverageMs',        1, ...
    'maxRs',                30, ...       % MOhm, on min(Rs, Rs_headerString)
    'fallbackMaxRs',        35);          % only when exactly one sweep passes maxRs

end
