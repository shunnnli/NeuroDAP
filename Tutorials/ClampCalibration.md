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

Calibration uses **all recorded LabJack channels** in `labjack.raw`, with
channel names and modulation settings from `labjack.name`, `labjack.mod`, and
`labjack.modFreq`. The NI `photometry_raw` signal is not used. Inputs are
`labjack`, `redClamp`, and `blueClamp` in `data_*.mat`, plus `params` in
`sync_*.mat`. Run `loadSessions` first to load and synchronize these recordings.
Behavioral trials, camera/licking fields, and `timeseries_*.mat` are not needed.

Laser events are detected on the NI recording. Fluorescence is interpolated
onto those NI samples using the full `params.sync.timeNI` and
`params.sync.timePhotometry` vectors produced by `loadSessions` /
`assignTimeStamp`. This handles different start times, sample rates, recording
lengths, and relative clock drift; it does not assume the devices started
together. Missing or invalid sync timestamps cause an error asking for
synchronized inputs. Data outside the LabJack recording are not extrapolated.

Modulated channels first use `demodulateSignal`'s `demodData_nodetrend` output;
its spectral-window center timestamps are mapped through the LabJack sync
clock. Unmodulated channels use their recorded fluorescence voltage. This
preserves fluorescence baselines instead of dividing already normalized
rolling-z or dF/F signals by a near-zero baseline. No Arduino ADC conversion or
quantization is applied to LabJack data.

Following `BrainClamp/PID-tuning.ipynb`, the aligned signal is sampled at
approximately 200 Hz and receives a 20 ms exponential filter. Each trial is
normalized as `(F-F0)/F0`, with `F0` the mean during the preceding second.
The saved analysis records the actual sampling rate and resolved options.
Existing MAT variables are preserved when analysis, events, and settings are
saved.

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
baseline windows, incomplete LabJack coverage during stimulation, and zero
baselines are excluded from response statistics; missing post-stimulation
samples are represented as NaN.

Each LabJack channel gets its own figure showing mean ± SEM photometry traces
for every power, separately for red and blue. Dotted lines mark each group's
median stimulation offset. The
lower panels show individual trial responses, their median at each power,
and a linear fit. Each trial response is the median ΔF/F0 over its last
2 seconds of stimulation (or the full pulse for shorter stimuli).

Outputs in the session directory:

- `analysis_*.mat`: `calibrationAnalysis.signals(k)` holds channel `k`'s name,
  alignment metadata, events, power groups, and fits. Group fields include
  `raw_signal`, `baseline_signal`, and `dff`; raw units are volts or demodulated
  amplitude. `calibrationAnalysis.signalNames` lists the channel order.
- `behavior_*.mat`: `calibrationEvents`; original NI indices are one-based,
  with offsets pointing to the first sample after stimulation. `onset_sync_sec`
  uses the shared sync clock; `onset_sec` is relative to NI recording start.
  `validTrial` has one column per LabJack channel in `signalNames` order.
- `Calibration_photometry_power_response_01_<channel>/` (and subsequent
  channel numbers): PNG and PDF plots, with the channel name in each title.

For example, `result.signals(1).groups` contains the first LabJack channel's
per-power traces and statistics. Each signal's `events.validTrial` records its
own usable trials; the NI laser events remain the same across signals.

Set `plotPhotometry=false` to save the analysis without plotting. Use
`redo=false, analyzeTraces=false` to reuse saved results when calibration
settings match. Legacy NI-only caches are automatically recomputed using
LabJack. Calibration bypasses behavioral trials and clamp-target processing;
no PID configuration files are edited.
