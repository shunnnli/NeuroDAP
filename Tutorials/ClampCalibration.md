# Clamp calibration sessions

`analyzeSessions_clamp` automatically selects calibration mode when the session
name contains `calibration` (case-insensitive), even if a behavioral task was
previously saved. It can also be selected with `task='calibration'`.
All calibration session handling is delegated to
`analyzeSessions_clampCalibration`, which can also be called directly and
returns the analysis structure. It runs in four numbered steps that match the
code blocks in the file:

1. **Load mat files** — `params` from `sync_*.mat`, the clamp commands from
   `data_*.mat`, the photometry channels from `timeseries_*.mat`.
2. **Detect clamp pulses** in the NI commands, with `findClampPulses`.
3. **Group pulses by laser power**, keeping `onsetIdx`, the starting NI sample
   of every pulse at that power.
4. **Align and plot** — those starting samples go into `plotTraces`, once per
   channel and power level, followed by the late-on response and its power fit.

```matlab
analyzeSessions_clamp(sessionpath);

% Run calibration directly and return the grouped responses.
result = analyzeSessions_clampCalibration(sessionpath);

% Override DAC ranges to match the settings used during acquisition.
analyzeSessions_clamp(sessionpath,blueClampRange=[800 1500],redClampRange=[25 350]);

% Optional calibration-specific settings (seconds, except powerTolerance,
% thresholdPct, targetFs, and baselineSubtract).
analyzeSessions_clamp(sessionpath,calibration=struct( ...
    'preTime',1,'postTime',2,'lateTime',2,'maxDuration',5,'powerTolerance',2));
```

Supported `calibration` fields are `targetFs`, `thresholdPct`, `minDuration`,
`maxDuration`, `minGap`, `powerTolerance`, `preTime`, `postTime`, `lateTime`,
and `baselineSubtract`. Unrecognized fields raise an error instead of being
silently ignored.

Calibration analyzes the **processed photometry in `timeseries_*.mat`**: every
`timeSeries` channel whose `system` is `NI` or `LJ`, using the stored `data`,
`finalFs`, and `system` fields. The `blueClamp`, `redClamp`, and `clampTarget`
command channels are skipped, as are camera channels — they are inputs to the
sweep, not responses. `labjack.raw` and `photometry_raw` are no longer read.
Inputs are `timeseries_*.mat`, `redClamp` and `blueClamp` in `data_*.mat`, and
`params` in `sync_*.mat`. Run `loadSessions` first to process and synchronize
these recordings. Behavioral trials, camera/licking fields, and `clampTarget`
are not needed.

Laser events are detected on the NI clamp commands, so pulse indices are NI
samples. `plotTraces`/`getTraces` then map each onset through
`params.sync.timeNI` onto the channel's own clock: the stored trace is treated
as uniformly sampled at `finalFs` starting at `params.sync.timePhotometry(1)`
for LabJack channels and `params.sync.timeNI(1)` for NI channels. This is the
same alignment every other session analysis uses, and it corrects start-time
offsets and different sampling rates. Unlike per-sample interpolation, it does
not track relative clock drift within a session. Missing `timeNI`, missing
`timePhotometry` for a LabJack channel, or a `timeNI` length that disagrees with
the NI recording cause an error asking for synchronized inputs.

Traces are plotted in the units `loadSessions` stored them in, recorded per
channel in `signalUnits` (`z-score` for `rolling-z`, `\DeltaF/F` for `dff`).
The mean of the pre-onset window is subtracted for the response statistics,
which can be turned off with `baselineSubtract=false`. No second normalization
is applied: already normalized rolling-z or dF/F signals are never divided by a
near-zero baseline. There is no exponential filter and no Arduino ADC conversion
or quantization — the stored trace is already downsampled and detrended.

The saved analysis records the detection settings, the resolved options, the
detection sampling rate, and each channel's `finalFs`. Existing MAT variables
are preserved when analysis, events, and settings are saved.

Red is labeled excitation and blue inhibition. Calibration DAC defaults are
`[25 500]` for red and `[800 1600]` for blue, as in `BrainClamp/PID-tuning.ipynb`;
other session types retain their existing defaults. Powers denote percent of
these DAC ranges, not measured optical power. Previously saved percentage traces
are not reused, so changing the ranges updates power estimates.

Pulses exceed 5% command, last 0.5–5 seconds, and merge off gaps of up to
2 seconds, following the notebook's detector, applied to the commands decimated
to approximately 200 Hz. Unlike its strict `<5 s` filter, exactly 5-second
pulses are retained. Powers within 2 percentage points are grouped and labeled
with their rounded median in `level_pct`. Distinct levels are retained even with
few repeats. Incomplete stimulation, missing full baseline windows, and
incomplete photometry coverage during stimulation are excluded from response
statistics; missing post-stimulation samples are represented as NaN.

Each photometry channel gets its own 2x2 figure. The top row is `plotTraces`
output — mean ± SEM per power, red on the left and blue on the right, with more
opaque lines for higher powers and a patch marking the stimulation window. The
bottom row shows individual pulse responses, their median at each power, and a
linear fit. Each pulse response is the median over its last 2 seconds of
stimulation (or the full pulse for shorter stimuli).

Outputs in the session directory:

- `analysis_*.mat`: `calibrationAnalysis.powerGroups` is the step-3 grouping
  shared by every channel — `channel`, `power_pct`, `nPulses`, `pulseRows`,
  `onsetIdx`, `onsetTime_sec`, `duration_sec`, and `complete`.
  `calibrationAnalysis.signals(k)` holds channel `k`'s name, system, units,
  detrending provenance, alignment metadata, fits, and its own copy of the
  groups with the step-4 results appended: `traces` (as `plotTraces` returned
  them), `baseline` (the per-pulse pre-onset mean), `response`
  (baseline-subtracted), `mean_signal`, `sem_signal`, `lateOnMedian`,
  `medianResponse`, `nTrials`, and `validPulse`.
  `calibrationAnalysis.signalNames` lists the channel order, and
  `calibrationAnalysis.detection` records the shared pulse detection.
- `behavior_*.mat`: `calibrationEvents`; original NI indices are one-based,
  with offsets pointing to the first sample after stimulation. `onset_sync_sec`
  uses the shared sync clock; `onset_sec` is relative to NI recording start.
  `validTrial` has one column per photometry channel in `signalNames` order.
- `Calibration_photometry_power_response_01_<channel>.pdf` (and subsequent
  channel numbers): vector PDF plots, with the channel name in each title.

For example, `result.signals(1).groups` contains the first photometry channel's
per-power traces and statistics. Each signal's `events.validTrial` records its
own usable trials; the NI laser events remain the same across signals.

Set `plotPhotometry=false` to save the analysis without plotting; the traces are
still extracted, just not drawn. `redo` and `analyzeTraces` are accepted so
`analyzeSessions_clamp` can forward them, but the calibration analysis always
recomputes — the photometry is already processed, so there is nothing expensive
to cache. Calibration bypasses behavioral trials and clamp-target processing; no
PID configuration files are edited.
