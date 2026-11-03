# DMD hotspot analysis: parity with the pyNeuroDAP notebook

**Date:** 2026-11-02
**Status:** Draft, awaiting approval

## Problem

`Shun_loadSliceData_DMD.m` and the notebook at
[pyNeuroDAP `notebooks/Shun_dmd_analysis.ipynb`](https://github.com/shunnnli/pyNeuroDAP/blob/main/notebooks/Shun_dmd_analysis.ipynb)
analyze the same recordings. The notebook reads the `.mat` files `loadSlicesDMD.m` writes, so there
is one dataset and two pipelines — but they disagree about which spots are hotspots, and the MATLAB
figures therefore highlight a different set of spots than the figures in the paper analysis.

The disagreement is not in the baseline levelling, which is already per sweep per spot on both
sides. It is in the hotspot test, in three separate ways:

1. **The offline re-test differs on five points** (`loadSlicesDMD.m:399-404` against the notebook's
   `test_trace_for_hotspot`): MATLAB judges one sweep at a time where the notebook judges the
   spot's sweep average; MATLAB's noise SD comes from the sweep-wide baseline window and is shared
   by every spot in the sweep where the notebook takes each spot's own pre-stimulus window; MATLAB
   tests 0–30 ms where the notebook tests 0–50 ms; MATLAB has no amplitude floor where the notebook
   requires a 2 pA peak; and MATLAB ORs the verdict across repetitions
   (`loadSlicesDMD.m:761`) where the notebook averages first.
2. **The notebook's operative definition is not that test at all.** Its live configuration sets
   `HOTSPOT_TYPE = "evoked"`, a charge-excess detection metric MATLAB has no equivalent of.
3. **`timeRange = [-10,50]`** in the current script truncates the stored per-spot window so far that
   the notebook's own windows do not fit: its −30 ms baseline clips to −10 ms, and
   `detection_windows` raises `no 45 ms control window fits after 50 ms`.

Separately, MATLAB applies no series-resistance filter, where the notebook drops sweeps above
30 MΩ before computing anything.

## Goals

- The spots drawn as hotspots in the MATLAB search figures are the spots the notebook analyzes.
- Their complement is the null-spot set used by the statistics panels in those figures.
- Spots the notebook would not analyze remain **visible**, drawn in a more transparent colour,
  rather than being dropped from the figure.
- One baseline window, −30 to −1 ms, serving both the levelling and the noise estimate.
- The hotspot determination lives with the rest of the analysis, in `loadSlicesDMD.m`.
  `analyzeDMDSearch.m` consumes a verdict; it does not compute one.

## Non-goals

- Changing the notebook or anything else in the pyNeuroDAP repository. None is needed — the
  notebook already uses the `(-30, -1)` ms baseline this design adopts.
- Rewriting `analyzeSlice_DMD.m`'s cross-cell aggregation, `analyzeDMDSearchPair.m`, or the
  `Other methods/ScriptsLuca/` forks.
- Reproducing the notebook's cross-cell statistics (contributions, concentration, group
  comparisons) in MATLAB. Only the hotspot verdict is brought over.
- Reproducing the notebook's **search selection** — epoch-name parsing, Vhold hints from name
  tokens, and the `DEPTH_CHOICE` marker files. See "Decisions made".
- Changing the stored per-sweep `responses.hotspot` flag. See "Decisions made".

## Approach

Everything the verdict needs already exists inside `loadSlicesDMD.m` at the point where a search's
depths are aggregated: `depthCurrentMap{spotIdx}` is a reps × samples matrix, `spotsAtDepth`
carries `Sweep`, `Repetition` and the per-sweep `QC` struct, and the cell boundary that the
cell-level control reference needs is already tracked there.

So the verdict is computed once, at load time, and stored per spot. `analyzeDMDSearch.m` reads it.
The arithmetic lives in one new function so it can be tested on its own and so the two verdicts
stay in one place.

Approaches rejected:

- **Computing the verdict in `analyzeDMDSearch.m` at plot time.** An earlier draft did this, on the
  grounds that it avoided touching stored data and kept the criteria switchable without a reload.
  But all three §2 corrections change stored data anyway, so the reload is not avoidable; it would
  have pushed the Rs join, the cell-reference pass and the sweep bookkeeping into a plotting
  function, and it would have needed a new stored field (`searchRepMap`) purely to recover the
  sweep identity that `spotsAtDepth` already has at load time. Storing the *metrics* rather than
  the boolean recovers most of the switchability — see below.
- **Porting the notebook's full selection pipeline into MATLAB.** The notebook's filters are
  cross-cell and cross-animal (`ANIMALS_TO_EXCLUDE`, group assignment, `DEPTH_CHOICE` markers);
  `analyzeDMDSearch` is a single-cell single-search function with no visibility of them.

### 1. `Shun_loadSliceData_DMD.m` — split storage and plotting windows

Today one `timeRange` variable (line 21) feeds `loadSlices`, `loadSlicesDMD` **and**
`analyzeSlice_DMD` (line 64), which passes it to `analyzeDMDSearch` where it sets each tile's
`xlim` and the control-data slicing. Storage and display are the same number.

They have to separate, because they want different values:

- **`storageTimeRange = [-30,100]`** for `loadSlicesDMD`. The baseline is −30 to −1 ms and the
  notebook can only see what the stored window contains, so the window must start at or before
  −30 ms. The +100 ms end is what makes the detection control window (55–100 ms) exist at all.
- **`plotTimeRange = [-10,50]`** for `analyzeSlice_DMD`, matching the notebook's
  `PLOT_WINDOW_MS = (-10.0, 50.0)`, which is display-only there and affects no metric.

Without the split, taking `[-30,100]` literally would silently re-scale every per-spot tile from
60 ms wide to 130 ms, with 30 ms of baseline and 50 ms of tail crowding out the response.
`analyzeSlice_DMD` needs a second option to carry the plotting value through.

### 2. `Methods/loadSlicesDMD.m` — the main change

**Two corrections to what is already there:**

- **Control window → −30 to −1 ms.** Currently
  `controlWindow = stimOnset - controlWindowSamples - 1 : stimOnset - 1`, which is −30.1 to −0.1 ms
  and runs right up to the pulse. Becomes
  `stimOnset - controlWindowSamples : stimOnset - guardSamples` with
  `guardSamples = round(options.baselineGuard * outputFs/1000)` and a new
  `options.baselineGuard double = 1` (ms). This is the 1 ms guard band the notebook uses; it keeps
  the leading edge of the stimulus artifact and any capacitive transient out of the baseline. Both
  the levelling median and `responses.control` follow from it.
- **`analysisWindow` off-by-one.** `eventSample = abs(outputFs*timeRange(1)/1000)` at
  `loadSlicesDMD.m:304` indexes one sample before 0 ms, since index 1 is `timeRange(1)`. Add 1.

**And the new work, in the search-aggregation block (`loadSlicesDMD.m:634-830`):**

After the depth loop has filled `depthCurrentMap` for every depth of a search, two passes:

1. **Build the search's control reference** from every spot's Rs-filtered sweep average.
2. **Score each spot** against it.

Per spot, this stores alongside the existing maps:

| Field | Contents |
| --- | --- |
| `searchChargeExcessZ{depthIdx}{spotIdx}` | charge excess, cell-null z |
| `searchChargeExcessPc{depthIdx}{spotIdx}` | the same as pC over the window |
| `searchDerivativeExcessZ{depthIdx}{spotIdx}` | path-length excess, z |
| `searchReanalysisHotspot{depthIdx}{spotIdx}` | the corrected offline re-test, per spot |
| `searchNSweepsKept{depthIdx}{spotIdx}` | sweeps surviving the Rs filter |

Storing the **metrics** rather than a boolean is what keeps the criterion switchable:
`DETECTION_CRITERION` and `minChargeZ` are applied at plot time, so changing either needs no
reload. Changing a *window* does need one, which is acceptable — those are not knobs that get
turned.

**Rs filter.** Rows of `depthCurrentMap{spotIdx}` whose sweep fails the limit are dropped before
the spot is averaged. The limit is the notebook's `QC_METRIC_CHOICE = "min"`: judge
`min(Rs, Rs_headerString)` against 30 MΩ, treating a non-positive `Rs_headerString` as missing
rather than as an unbeatably small Rs. Both halves matter — on this dataset, 35-on-`Rs` and
30-on-`min` keep substantially different sweep sets, enough to make a parity comparison fail for
reasons unrelated to the port. The join is direct: `spotsAtDepth` already carries `Sweep`,
`Repetition` and `QC` per row, so no new identity field is needed.

**Rs-filtered statistics.** A parallel `cellStats.filtered.*` holding `max`, `min`, `maxTime`,
`minTime` and `auc` computed the notebook's way — the extremes of the **Rs-filtered sweep mean**,
each averaged over `PEAK_AVERAGE_MS` = 1 ms around its own sample. The existing `cellStats.*`
(per-sweep peak-window averages, then averaged across reps, over every sweep) is left untouched so
`analyzeSlice_DMD` and `analyzeDMDSearchPair` keep working. The two are different quantities, not
two estimates of one.

The stale comment at `loadSlicesDMD.m:297` ("50ms before each spot stim onset") is corrected to
match `controlWindowLength`, which is 30.

**Every result folder must be regenerated.** All of the above changes stored data, and existing
`Results-*` folders only hold a 20 ms pre-stimulus window, so the −30 ms baseline cannot be taken
from them at all. `analyzeDMDSearch` requires the new fields and errors with a message naming the
re-run, rather than carrying a fallback path for old files.

### 3. `Methods/dmdHotspotVerdict.m` — new

One function, both notebook definitions, no plotting, no file access. Two entry points, because the
detection metrics need a cell-level reference that one spot cannot supply:

```
ref     = dmdHotspotVerdict('reference', spotTraces, timeMs, criteria)
verdict = dmdHotspotVerdict('score',     spotTraces, timeMs, criteria, ref)
```

**Reanalysis verdict** — the notebook's `test_trace_for_hotspot` with `trace_mode = "mean"`:
average the reps; take the SD (normalized by N−1) of the average over the baseline window; a spot
passes when at least `minSamples` samples of `|average|` inside `windowMs` reach
`thresholdFactor × SD`, and the window's peak `|current|` reaches `minAmplitudePa`.

**Detection metrics** — the notebook's `spot_detection_metrics`:

- Response window 5–50 ms. Control windows are consecutive non-overlapping windows of the same
  length starting 5 ms after it ends, i.e. 55–100 ms.
- `chargeDensity(W)` = sum of `|trace − median(trace over W)|` over W, × sample interval in ms
  / 1000, / duration in ms. The median subtraction ("de-levelling") is per window.
- `derivativeDensity(W)` = sum of `|trace(t+lag) − trace(t)|` over W / duration, `lag` = 2 ms.
- The reference is the median of each density over the control windows, taken across the spots of
  **one depth**; the spread is `1.4826 × MAD` of that same set, `NaN` below 5 spots.
- `chargeExcessZ = (chargeDensity(response) − reference.chargeCentre) / reference.chargeSpread`,
  likewise for the derivative.
- A spot is an evoked hotspot when `chargeExcessZ >= minChargeZ`.

Both run on the stored trace with each sweep's own baseline median removed
(`baseline_corrected_sweeps`), which on top of the corrected levelling is close to a no-op but
keeps the two pipelines arithmetically identical.

Defaults mirror the notebook's live configuration:

| Field | Value |
| --- | --- |
| `baselineMs` | `[-30, -1]` |
| `thresholdFactor` | 3 |
| `minSamples` | 50 |
| `reanalysisWindowMs` | `[0, 50]` |
| `minAmplitudePa` | 2 |
| `detectionWindowMs` | `[5, 50]` |
| `detectionGapMs` | 5 |
| `derivativeLagMs` | 2 |
| `minChargeZ` | 0.5 |
| `minSpotsForReference` | 5 |
| `maxRs` | 30 MΩ, on `min(Rs, Rs_headerString)` |

### 4. `Methods/analyzeDMDSearch.m` — consume the verdict

Deliberately small. A new `options.hotspotType string = "evoked"`, accepting `"acquisition"`,
`"reanalysis"` and `"evoked"`, plus `options.minChargeZ` and `options.nullAlpha`.

- **Verdict.** `hotspotSpot = chargeExcessZ >= options.minChargeZ` for `"evoked"`;
  `searchReanalysisHotspot` for `"reanalysis"`. This replaces the `sum(x)>=1` OR at
  `analyzeDMDSearch.m:273-274`, where the existing `% (Optional) Change hotspot criteria` comment
  already anticipates it.
  Note that OR currently reads `responseMap.hotspot`, the computed per-sweep test — *not*
  `isResponse`. `"acquisition"` mode must therefore source the DMD software's own call from
  `responseMap.isResponseMap`, by looking up each spot's box through `search_spotLocation` the way
  the reconstruction path already does (`loadSlicesDMD.m:2163`).
- **Transparency.** Spot patches and traces for hotspots keep their current colour at full opacity.
  Everything else — spots failing the verdict, and spots whose sweeps were all dropped by the Rs
  filter — is drawn in the same colour at `options.nullAlpha` (default 0.3). Nothing is removed
  from the figure.
- **Statistics panels.** The hotspot/null split feeding the min/max, AUC, absolute AUC, latency and
  E/I rate panels (`analyzeDMDSearch.m:408-508`) uses the verdict, and the values come from
  `cellStats.filtered.*` so the bars and the traces above them agree about which sweeps are in.
- **Control-data slicing.** The new control window is 291 samples, down from 302, so the
  `ctrlData(:, end-length(analysisWindow)+1:end)` slice at `analyzeDMDSearch.m:252` would compute a
  negative start index and throw. The commented-out robust version at `:255-257` replaces it.

## Decisions made

These are the builder's calls, not the user's. Flagged so they can be checked while reading.

- **No scope subtitle.** An earlier draft proposed annotating whether the notebook would analyze
  this search at all. The real rule is `find_repeated_final_depth_searches`: the epoch name must
  contain `hotspot`, `maxSearchDepth` decides between two candidates, a Vhold hint parsed from
  `m35`/`exci`/`inhi` tokens overrides the measured value when they disagree by more than 15 mV,
  and a per-cell `DEPTH_CHOICE-MAX|SHALLOW|BOTH|NONE` marker file in the `cellN/` folder selects
  which depths count (`slice.py:224-253`, default `"max"`) via a `Depth choice` column that does
  not exist in `cells_DMD_*.mat`. A loose approximation is wrong in practice — CM023 cell 4 has
  three searches at `Vhold ≈ −34.9` with depths 3, 2, 3, of which a Vhold-and-depth rule marks two
  in scope where the notebook keeps one. Per-spot transparency already carries the answer that
  matters, so no scope claim is made.
- **The stored `responses.hotspot` is left alone.** The notebook never reads it under any
  `HOTSPOT_TYPE`: `"acquisition"` reads `isResponse`, and the other two recompute from traces. It
  never reads `nullSpotData` either. So rewriting it buys no parity, and the per-sweep flag stays
  as the legacy quantity it is. The noise panels keep that old definition and are labelled as such.

  This is *not* an argument that the noise model stays fixed: `preStimData`
  (`loadSlicesDMD.m:420`) is accumulated from `controlWindow`, so the §2 control-window change
  moves the Gaussian behind `Ethres`/`Ithres` anyway. The dashed threshold lines and `getYlimit`
  y-ranges will shift slightly on every tile. That is accepted, not guarded against.
- **The derivative excess is computed but does not gate.** The notebook computes both and
  `DETECTION_CRITERION = "charge"` selects charge alone. Computing both costs nothing and makes the
  criterion switchable without rescoring.
- **The control reference is per depth, falling back to the search.** Each spot contributes one
  reading — the median charge density over its 55–100 ms control windows — and the median and
  scaled MAD of that depth's readings are what every spot at that depth is scored against. When a
  depth holds fewer than `minSpotsForReference` = 5 spots, the reference is rebuilt from every spot
  of the search instead; if the search is also under 5, the spread is `NaN`, the z is `NaN`, and no
  spot there is a hotspot.

  The fallback cannot cost parity, which is why it is worth having. The notebook only ever scores
  the repeated final-depth search, and that search has far more than 5 spots, so the fallback does
  not fire on any spot the notebook scores. In the one case where it could — a repeated search with
  under 5 spots — that search holds a single depth, so the search-wide pool is the same set and the
  result is `NaN` either way. Every depth the fallback actually rescues is a coarse search depth
  the notebook never looks at, where a scored verdict is more useful on the figure than a blank one.

  Strictly, the notebook keys the reference on `CELL_KEYS = [animal, results_folder, cell]` —
  depth is not in the key. It coincides with per-depth because
  `find_repeated_final_depth_searches` returns one repeated search per cell under the default
  `DEPTH_CHOICE-MAX`, and `repeated_search_depth` gives that search exactly one depth. Under
  `DEPTH_CHOICE-BOTH`, or if a cell has two qualifying `maxSearchDepth` searches at −35 mV, the
  notebook pools them under one cell key and this design would not. Nothing in the dataset is
  known to hit that case; if the parity test shows a systematic z offset on one cell, this is
  where to look.
- **`cellStats.filtered.*` is additive.** The existing `cellStats.*` is left in place rather than
  overwritten, so `analyzeSlice_DMD` and `analyzeDMDSearchPair` are untouched by this work.

## Resolved: inter-spot interval

An earlier draft flagged that the stored window running to +100 ms would overlap the next spot's
stimulus if the inter-spot interval were under 100 ms, which would put a real response inside the
55–100 ms control windows and bias every charge excess toward calling real hotspots null.

`protocol.isi` is 200 ms uniformly across the searches checked (CM023, SL458, SL468, CM014), so the
control windows sit well inside the gap. No mitigation is needed. Worth a one-line assertion in the
code so a future dataset with a shorter ISI fails loudly rather than silently biasing.

## Verification

This is the acceptance test, and it can be run today — nothing is blocked on another repository:

- A search re-run end to end, with the per-spot `chargeExcessZ` from MATLAB compared against the
  notebook's `charge_excess_z` for the same cell, search and depth. They should agree to numerical
  precision. This is the test that the port is faithful.
- The hotspot set from `hotspotType = "evoked"` compared against the notebook's `is_hotspot` on the
  same final-depth search.
- `hotspotType = "acquisition"` reproduces today's `isResponseMap` spot for spot, confirming the
  pixel-map lookup is right. This checks the verdict only — the bar values change under every mode,
  since they now come from `cellStats.filtered.*`.
- A depth-1 search still renders, with every spot transparent and no empty-axis errors from the
  statistics panels.
- A sweep with `Rs_headerString <= 0` is judged on its computed `Rs` alone, not dropped and not
  treated as passing.
- A spot whose sweeps are all dropped by the Rs filter renders transparent and contributes to
  neither the control reference nor the statistics panels.
