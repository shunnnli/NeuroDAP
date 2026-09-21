# Clamp calibration sessions

`analyzeSessions_clamp` automatically selects calibration mode when the session
name contains `calibration` (case-insensitive), even if a behavioral task was
previously saved. It can also be selected with `task='calibration'`.
All calibration session handling is delegated to
`analyzeSessions_clampCalibration`, which can also be called directly and
returns the analysis structure. It uses `analyzeClampCalibration` for response
extraction and `plotClampCalibration` for plotting.

```matlab
analyzeSessions_clamp(sessionpath);

% Run calibration directly and return the grouped responses.
result = analyzeSessions_clampCalibration(sessionpath);

% Override DAC ranges to match the settings used during acquisition.
analyzeSessions_clamp(sessionpath,blueClampRange=[800 1500],redClampRange=[25 350]);

% Optional calibration-specific settings (seconds, except powerTolerance).
analyzeSessions_clamp(sessionpath,calibration=struct( ...
    'preTime',1,'postTime',2,'lateTime',2,'maxDuration',5,'powerTolerance',2));
```

Calibration uses `photometry_raw`, `redClamp`, and `blueClamp` from the NI
recording, at `params.sync.behaviorFs`. Following `BrainClamp/PID-tuning.ipynb`,
it samples at approximately 200 Hz, converts photometry voltage to rounded
10-bit ADC units, and applies a 20 ms exponential filter. Each trial is
normalized as `(F-F0)/F0`, with `F0` the mean during the preceding second.
The saved analysis records the actual sampling rate and all resolved options.
Only `data_*.mat` and `sync_*.mat` are required; the calibration path does not
load behavioral timeseries or require camera/licking fields. Existing MAT
variables are preserved when analysis, events, and settings are saved.

Red is labeled excitation and blue inhibition. Calibration DAC defaults are
`[25 500]` for red and `[800 1600]` for blue, as in the notebook; other session
types retain their existing defaults. Powers denote percent of these DAC
ranges, not measured optical power. Previously saved percentage traces are
not reused, so changing the ranges updates power estimates.

Pulses exceed 5% command, last 0.5–5 seconds, and merge off gaps of up to
2 seconds, following the notebook's detector. Unlike its strict `<5 s`
filter, exactly 5-second pulses are retained. Powers within 2 percentage
points are grouped and labeled with their rounded median. Distinct levels
are retained even with few repeats. Incomplete stimulation, missing full
baseline windows, and zero baselines are excluded from response statistics;
missing post-stimulation samples are represented as NaN.

The figure shows mean ± SEM photometry traces for every power, separately for
red and blue. Dotted lines mark each group's median stimulation offset. The
lower panels show individual trial responses, their median at each power,
and a linear fit. Each trial response is the median ΔF/F0 over its last
2 seconds of stimulation (or the full pulse for shorter stimuli).

Outputs in the session directory:

- `analysis_*.mat`: `calibrationAnalysis`, including detected events, aligned
  raw/ΔF/F0 trials, per-power summaries, fit slopes/intercepts/R², and settings.
- `behavior_*.mat`: `calibrationEvents`; original NI indices are one-based,
  with offsets pointing to the first sample after stimulation.
- `Calibration_photometry_power_response/`: PNG and PDF plots.

Set `plotPhotometry=false` to save the analysis without plotting. Use
`redo=false, analyzeTraces=false` to reuse saved results when calibration
settings match. Calibration bypasses behavioral trials and clamp-target
processing; no PID configuration files are edited.
