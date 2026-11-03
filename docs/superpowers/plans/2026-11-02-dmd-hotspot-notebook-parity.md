# DMD hotspot notebook parity — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the hotspots drawn opaque in the MATLAB DMD search figures be exactly the spots the pyNeuroDAP notebook analyzes, with their complement as the null set.

**Architecture:** Five small pure functions carry the notebook's arithmetic and are unit-tested on synthetic traces. `loadSlicesDMD.m` calls them in its search-aggregation block and stores per-spot metrics (not booleans, so thresholds stay switchable). `analyzeDMDSearch.m` reads those metrics, applies the threshold, and feeds the existing opacity switch and statistics split.

**Tech Stack:** MATLAB R2024b. Tests use `functiontests(localfunctions)` with `verifyEqual`/`verifyTrue`, matching `Methods/tests/testDmdProtocolClassification.m`.

MATLAB is **not on `PATH`** on this machine — it lives at `/Applications/MATLAB_R2024b.app/bin/matlab`, verified working in `-batch` mode. Every `matlab -batch ...` command below is shorthand for that full path. Run from the repository root so `Methods/` resolves; prepend `addpath(genpath('Methods'));` inside the `-batch` string if a test cannot find a function.

**Spec:** [docs/superpowers/specs/2026-11-02-dmd-hotspot-notebook-parity-design.md](../specs/2026-11-02-dmd-hotspot-notebook-parity-design.md)

## Global Constraints

- Baseline window: `[-30, -1]` ms, used for both levelling and the noise SD.
- Reanalysis test: `thresholdFactor` 3, `minSamples` 50, window `[0, 50]` ms, `minAmplitudePa` 2, on the sweep mean, SD normalized by N−1.
- Detection: window `[5, 50]` ms, `detectionGapMs` 5, `derivativeLagMs` 2 ms, `minChargeZ` 0.5, `minSpotsForReference` 5.
- Rs limit: `min(Rs, Rs_headerString) <= 30` MΩ; a non-positive `Rs_headerString` is missing, not a small Rs.
- `peakAverageMs` 1 ms, centred on the extremum's own sample.
- Storage window `[-30,100]` ms; plot window `[-10,50]` ms.
- New files go in `Methods/`, tests in `Methods/tests/`, following the existing `arguments` block style.
- All densities are per ms of window; charge in pC, path length in pA.

## Review Focus

Input classes the spec implies but no task's happy-path test exercises. Each has a test added to the task that owns the code.

1. **A spot whose sweeps are all dropped by the Rs filter** — empty trace matrix must yield `NaN` metrics and render transparent, not error. (Task 6, Task 8)
2. **`Rs_headerString <= 0`** — judged on computed `Rs` alone; neither dropped nor auto-passed. (Task 3)
3. **A depth with fewer than 5 spots** — reference falls back to the search; if that is also under 5, spread is `NaN`, z is `NaN`, no hotspot. (Task 5)
4. **A trace too short for a control window** — must raise a named error rather than score against a truncated window. (Task 1)
5. **`NaN` samples inside a window** — ignored by the density and by the baseline SD, not propagated into the whole metric. (Task 2, Task 4)

---

### Task 1: Detection windows

**Files:**
- Create: `Methods/dmdDetectionWindows.m`
- Test: `Methods/tests/testDmdDetectionWindows.m`

**Interfaces:**
- Produces: `[windowMs, controlWindows] = dmdDetectionWindows(timeMs, windowMs, gapMs)` — `windowMs` a 1×2 `[start end]`, `controlWindows` an N×2 matrix of consecutive non-overlapping windows of the same duration, the first starting at `windowMs(2) + gapMs`.

- [ ] **Step 1: Write the failing tests**

```matlab
function testControlWindowsTileAfterGap(testCase)
timeMs = linspace(-30,100,1301);
[w,c] = dmdDetectionWindows(timeMs,[5 50],5);
verifyEqual(testCase,w,[5 50]);
verifyEqual(testCase,c,[55 100],'AbsTol',1e-9);
end

function testSecondWindowAppearsWhenRoomAllows(testCase)
timeMs = linspace(-30,200,2301);
[~,c] = dmdDetectionWindows(timeMs,[5 50],5);
verifyEqual(testCase,size(c,1),3);
verifyEqual(testCase,c(2,:),[100 145],'AbsTol',1e-9);
end

function testTooShortTraceErrors(testCase)   % Review Focus 4
timeMs = linspace(-10,50,601);
verifyError(testCase,@() dmdDetectionWindows(timeMs,[5 50],5), ...
            'dmdDetectionWindows:noControlWindow');
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdDetectionWindows.m')"`
Expected: FAIL, undefined function `dmdDetectionWindows`.

- [ ] **Step 3: Implement `dmdDetectionWindows`**

Walk a cursor from `windowMs(2)+gapMs`, emitting a window while `cursor + duration <= max(timeMs)` (with a small tolerance). Error id `dmdDetectionWindows:noControlWindow` when none fits.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdDetectionWindows.m')"`
Expected: PASS, 3 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdDetectionWindows.m Methods/tests/testDmdDetectionWindows.m
git commit -m "feat: add DMD detection window helper"
```

---

### Task 2: Window densities

**Files:**
- Create: `Methods/dmdWindowDensity.m`
- Test: `Methods/tests/testDmdWindowDensity.m`

**Interfaces:**
- Produces: `[chargePcPerMs, travelPaPerMs] = dmdWindowDensity(trace, timeMs, windowMs, derivativeLagMs)`. Charge de-levels by the window's own median before rectifying. Both ignore `NaN`. `trace` is one row vector.

- [ ] **Step 1: Write the failing tests**

```matlab
function testConstantOffsetContributesNothing(testCase)
timeMs = 0:0.1:45; trace = 7*ones(size(timeMs));
[q,g] = dmdWindowDensity(trace,timeMs,[0 45],2);
verifyEqual(testCase,q,0,'AbsTol',1e-12);
verifyEqual(testCase,g,0,'AbsTol',1e-12);
end

function testChargeDensityIsPerMillisecond(testCase)
timeMs = 0:0.1:45;
trace = zeros(size(timeMs)); trace(timeMs>=10 & timeMs<20) = -20;
[q,~] = dmdWindowDensity(trace,timeMs,[0 45],2);
verifyGreaterThan(testCase,q,0);          % rectified, so positive
end

function testNaNSamplesAreIgnored(testCase)   % Review Focus 5
timeMs = 0:0.1:45; trace = zeros(size(timeMs));
trace(timeMs>=10 & timeMs<20) = -20; clean = dmdWindowDensity(trace,timeMs,[0 45],2);
trace(5) = NaN;
verifyEqual(testCase,dmdWindowDensity(trace,timeMs,[0 45],2),clean,'RelTol',0.01);
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdWindowDensity.m')"`
Expected: FAIL, undefined function.

- [ ] **Step 3: Implement `dmdWindowDensity`**

Mask `timeMs` inclusively. Charge: subtract `median(...,'omitnan')`, `sum(abs(...),'omitnan')`, × sample interval in ms / 1000, ÷ duration in ms. Travel: `sum(abs(v(1+lag:end)-v(1:end-lag)),'omitnan')` ÷ duration, with `lag = max(1,round(derivativeLagMs/sampleIntervalMs))`.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdWindowDensity.m')"`
Expected: PASS, 3 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdWindowDensity.m Methods/tests/testDmdWindowDensity.m
git commit -m "feat: add DMD window charge and path-length densities"
```

---

### Task 3: Rs sweep filter

**Files:**
- Create: `Methods/dmdFilterSweepsByRs.m`
- Test: `Methods/tests/testDmdFilterSweepsByRs.m`

**Interfaces:**
- Produces: `keep = dmdFilterSweepsByRs(rs, rsHeader, maxRs)` — logical vector, same size as `rs`.

- [ ] **Step 1: Write the failing tests**

```matlab
function testJudgesTheMinimumOfBothEstimates(testCase)
keep = dmdFilterSweepsByRs([50 20 28],[25 40 29],30);
verifyEqual(testCase,keep,[true true true]);
end

function testDropsWhenBothEstimatesExceedLimit(testCase)
verifyEqual(testCase,dmdFilterSweepsByRs(50,45,30),false);
end

function testNonPositiveHeaderIsMissing(testCase)   % Review Focus 2
verifyEqual(testCase,dmdFilterSweepsByRs(45,0,30),false);
verifyEqual(testCase,dmdFilterSweepsByRs(20,-3,30),true);
end

function testNaNComputedFallsBackToHeader(testCase)
verifyEqual(testCase,dmdFilterSweepsByRs(NaN,25,30),true);
verifyEqual(testCase,dmdFilterSweepsByRs(NaN,NaN,30),false);
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdFilterSweepsByRs.m')"`
Expected: FAIL, undefined function.

- [ ] **Step 3: Implement `dmdFilterSweepsByRs`**

Mask `rsHeader <= 0` to `NaN`, take `min(rs, rsHeader)` ignoring `NaN`, compare `<= maxRs`. A sweep with no usable estimate does not pass.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdFilterSweepsByRs.m')"`
Expected: PASS, 4 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdFilterSweepsByRs.m Methods/tests/testDmdFilterSweepsByRs.m
git commit -m "feat: add Rs sweep filter matching notebook QC_METRIC_CHOICE"
```

---

### Task 4: Reanalysis hotspot test

**Files:**
- Create: `Methods/dmdDefaultCriteria.m`, `Methods/dmdHotspotVerdict.m`
- Test: `Methods/tests/testDmdHotspotVerdict.m`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `criteria = dmdDefaultCriteria()` — a struct carrying every value in Global Constraints, as the single source for both `loadSlicesDMD` and the tests. Fields: `baselineMs`, `thresholdFactor`, `minSamples`, `reanalysisWindowMs`, `minAmplitudePa`, `detectionWindowMs`, `detectionGapMs`, `derivativeLagMs`, `minChargeZ`, `minSpotsForReference`, `peakAverageMs`, `maxRs`.
  - `out = dmdHotspotVerdict('reanalysis', traces, timeMs, criteria)` where `traces` is reps × samples. Returns a struct with `baselineSd`, `threshold`, `nOverThreshold`, `peakAbs`, `isHotspot`.

- [ ] **Step 1: Write the failing tests**

```matlab
function testFlatTraceIsNotHotspot(testCase)
[timeMs,traces] = localSyntheticSpot(0);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());
verifyFalse(testCase,out.isHotspot);
end

function testLargeSustainedResponseIsHotspot(testCase)
[timeMs,traces] = localSyntheticSpot(-40);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());
verifyTrue(testCase,out.isHotspot);
end

function testAmplitudeFloorRejectsTinyButSignificant(testCase)
[timeMs,traces] = localSyntheticSpot(-1);   % clears 3*SD, under 2 pA
c = localCriteria(); c.minAmplitudePa = 2;
verifyFalse(testCase,dmdHotspotVerdict('reanalysis',traces,timeMs,c).isHotspot);
c.minAmplitudePa = 0;
verifyTrue(testCase,dmdHotspotVerdict('reanalysis',traces,timeMs,c).isHotspot);
end

function testSdComesFromBaselineWindowOnly(testCase)
[timeMs,traces] = localSyntheticSpot(-40);
out = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria());
verifyLessThan(testCase,out.baselineSd,5);   % response must not inflate it
end

function testNaNBaselineSamplesAreIgnored(testCase)   % Review Focus 5
[timeMs,traces] = localSyntheticSpot(-40);
clean = dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria()).baselineSd;
traces(1,10) = NaN;
verifyEqual(testCase, ...
    dmdHotspotVerdict('reanalysis',traces,timeMs,localCriteria()).baselineSd, ...
    clean,'RelTol',0.05);
end

function testDefaultCriteriaCarrySpecValues(testCase)
c = dmdDefaultCriteria();
verifyEqual(testCase,c.baselineMs,[-30 -1]);
verifyEqual(testCase,c.reanalysisWindowMs,[0 50]);
verifyEqual(testCase,c.minAmplitudePa,2);
verifyEqual(testCase,c.maxRs,30);
verifyEqual(testCase,c.minChargeZ,0.5);
end
```

Add local helpers `localCriteria()`, which returns `dmdDefaultCriteria()` so the tests and
`loadSlicesDMD` cannot drift apart, and `localSyntheticSpot(amplitudePa)` building 5 sweeps on
`timeMs = -30:0.1:100` with Gaussian noise (fixed `rng(0)`) plus a step of `amplitudePa` over
5–25 ms.

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdHotspotVerdict.m')"`
Expected: FAIL, undefined function.

- [ ] **Step 3: Implement `dmdDefaultCriteria` and the `'reanalysis'` branch of `dmdHotspotVerdict`**

`dmdDefaultCriteria` is a struct literal of the Global Constraints values. The `'reanalysis'`
branch means the sweeps (`'omitnan'`), takes `std(...,0,'omitnan')` over `criteria.baselineMs`,
counts samples of `abs(mean)` over `criteria.reanalysisWindowMs` reaching
`criteria.thresholdFactor * sd`, and requires both `>= criteria.minSamples` and
`peakAbs >= criteria.minAmplitudePa`.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdHotspotVerdict.m')"`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdDefaultCriteria.m Methods/dmdHotspotVerdict.m Methods/tests/testDmdHotspotVerdict.m
git commit -m "feat: add notebook reanalysis hotspot test and shared criteria"
```

---

### Task 5: Detection reference and scoring

**Files:**
- Modify: `Methods/dmdHotspotVerdict.m`
- Test: `Methods/tests/testDmdHotspotVerdict.m`

**Interfaces:**
- Consumes: `dmdDetectionWindows`, `dmdWindowDensity` from Tasks 1–2.
- Produces:
  - `ref = dmdHotspotVerdict('reference', spotTraces, timeMs, criteria)` where `spotTraces` is a cell array, one reps × samples matrix per spot. Returns `chargeCentre`, `chargeSpread`, `derivativeCentre`, `derivativeSpread`.
  - `out = dmdHotspotVerdict('score', traces, timeMs, criteria, ref)` returning `chargeExcessPc`, `chargeExcessZ`, `derivativeExcessPa`, `derivativeExcessZ`, `nSweeps`.

- [ ] **Step 1: Write the failing tests**

```matlab
function testResponsiveSpotScoresAboveFlatOnes(testCase)
timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);
flat = dmdHotspotVerdict('score',spots{1},timeMs,c,ref);
loud = dmdHotspotVerdict('score',localRespondingSpot(timeMs,-40),timeMs,c,ref);
verifyGreaterThan(testCase,loud.chargeExcessZ,flat.chargeExcessZ);
verifyGreaterThan(testCase,loud.chargeExcessZ,c.minChargeZ);
end

function testSpreadIsNaNBelowMinimumSpots(testCase)   % Review Focus 3
timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:3,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);
verifyTrue(testCase,isnan(ref.chargeSpread));
verifyTrue(testCase,isnan(dmdHotspotVerdict('score',spots{1},timeMs,c,ref).chargeExcessZ));
end

function testEmptyTracesScoreNaN(testCase)   % Review Focus 1
timeMs = -30:0.1:100; c = localCriteria();
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);
out = dmdHotspotVerdict('score',zeros(0,numel(timeMs)),timeMs,c,ref);
verifyTrue(testCase,isnan(out.chargeExcessZ));
verifyEqual(testCase,out.nSweeps,0);
end

function testSweepsAreBaselineCorrectedBeforeAveraging(testCase)
timeMs = -30:0.1:100; c = localCriteria();
base = localRespondingSpot(timeMs,-40);
offset = base + (1:size(base,1))'*25;   % per-sweep DC offsets
spots = arrayfun(@(k) localNoiseSpot(timeMs,k),1:8,'UniformOutput',false);
ref = dmdHotspotVerdict('reference',spots,timeMs,c);
verifyEqual(testCase, ...
    dmdHotspotVerdict('score',offset,timeMs,c,ref).chargeExcessZ, ...
    dmdHotspotVerdict('score',base,timeMs,c,ref).chargeExcessZ,'RelTol',1e-6);
end
```

Add `localNoiseSpot(timeMs,seed)` and `localRespondingSpot(timeMs,amplitudePa)` helpers.

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdHotspotVerdict.m')"`
Expected: FAIL on the four new tests; Task 4's six still pass.

- [ ] **Step 3: Implement the `'reference'` and `'score'` branches**

Both start by subtracting each sweep's own `median` over `criteria.baselineMs` and then averaging. `'reference'` takes each spot's median density over its control windows, then the median and `1.4826*mad(...,1)` across spots, with spread `NaN` below `criteria.minSpotsForReference`. `'score'` is `(responseDensity - centre)` for the excess, `÷ spread` for the z, `NaN` when the spread is not finite and positive.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdHotspotVerdict.m')"`
Expected: PASS, 10 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdHotspotVerdict.m Methods/tests/testDmdHotspotVerdict.m
git commit -m "feat: add charge-excess detection metrics and cell reference"
```

---

### Task 6: Sweep-mean statistics

**Files:**
- Create: `Methods/dmdSweepMeanStats.m`
- Test: `Methods/tests/testDmdSweepMeanStats.m`

**Interfaces:**
- Produces: `stats = dmdSweepMeanStats(traces, timeMs, windowMs, peakAverageMs)` with fields `max`, `min`, `maxTime`, `minTime`, `auc` — extremes of the baseline-corrected sweep **mean**, each averaged over `peakAverageMs` centred on its own sample; `auc` is `sum(mean)/Fs` over the window.

- [ ] **Step 1: Write the failing tests**

```matlab
function testPeakIsAveragedAroundItsOwnSample(testCase)
timeMs = -30:0.1:100; trace = zeros(size(timeMs));
trace(timeMs>=10 & timeMs<=12) = -30; trace(abs(timeMs-11)<0.05) = -100;  % one spike
s = dmdSweepMeanStats(trace,timeMs,[0 50],1);
verifyGreaterThan(testCase,s.min,-60);   % 1 ms average blunts the single sample
verifyEqual(testCase,s.minTime,11,'AbsTol',0.6);
end

function testEmptyTracesGiveNaN(testCase)   % Review Focus 1
s = dmdSweepMeanStats(zeros(0,1301),-30:0.1:100,[0 50],1);
verifyTrue(testCase,isnan(s.min) && isnan(s.max) && isnan(s.auc));
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "runtests('Methods/tests/testDmdSweepMeanStats.m')"`
Expected: FAIL, undefined function.

- [ ] **Step 3: Implement `dmdSweepMeanStats`**

Reuse the baseline-correct-then-average step from Task 5 (extract it to a shared local if convenient). Half-width in samples is `round(peakAverageMs/2/sampleIntervalMs)`, clamped to the window.

- [ ] **Step 4: Run to verify they pass**

Run: `matlab -batch "runtests('Methods/tests/testDmdSweepMeanStats.m')"`
Expected: PASS, 2 tests.

- [ ] **Step 5: Commit**

```bash
git add Methods/dmdSweepMeanStats.m Methods/tests/testDmdSweepMeanStats.m
git commit -m "feat: add sweep-mean peak statistics matching notebook definition"
```

---

### Task 7: `loadSlicesDMD` window corrections

**Files:**
- Modify: `Methods/loadSlicesDMD.m:19` (options), `:297` (comment), `:304-305` (analysis window), `:369` (control window)

**Interfaces:**
- Consumes: `dmdDefaultCriteria` from Task 4.
- Produces: `options.baselineGuard double = 1` (ms), `options.maxRs double = 30` (MΩ); `options.timeRange` default becomes `[-30,100]`. Any criteria not exposed as an option come from `dmdDefaultCriteria()`, so the Global Constraints live in one place.

- [ ] **Step 1: Change the control window to −30 to −1 ms**

At `:369`, replace `stimOnset - controlWindowSamples - 1 : stimOnset - 1` with
`stimOnset - controlWindowSamples : stimOnset - guardSamples`, where
`guardSamples = round(options.baselineGuard * options.outputFs/1000)` is computed beside
`controlWindowSamples` at `:300`. Add `options.baselineGuard double = 1` and
`options.maxRs double = 30` to the arguments block, change `options.timeRange` to `[-30,100]`, and
correct the stale `:297` comment from "50ms" to match `controlWindowLength`.

- [ ] **Step 2: Fix the analysis-window off-by-one**

At `:304`, `eventSample = abs(options.outputFs*options.timeRange(1)/1000) + 1`, so
`analysisWindow` starts at 0 ms rather than one sample before it.

- [ ] **Step 3: Assert the inter-spot interval leaves room for a control window**

Where `protocol.stimOnset` is available, error if the spacing is under
`detectionWindowMs(2) + detectionGapMs + diff(detectionWindowMs)` ms, naming `protocol.isi`, so a
future dataset with a shorter ISI fails loudly rather than scoring against the next spot's response.

- [ ] **Step 4: Verify the window arithmetic and that nothing else broke**

Run: `matlab -batch "runtests('Methods/tests')"`
Expected: PASS, including the existing suites. Then confirm the new control window is 291 samples
at 10 kHz and that `analysisWindow` now starts exactly at 0 ms:
`matlab -batch "assert(numel((3000-300):(3000-10))==291); assert(abs(((300+1)-301)/10)==0); disp('ok')"`
Expected: `ok`.

- [ ] **Step 5: Commit**

```bash
git add Methods/loadSlicesDMD.m
git commit -m "fix: align DMD baseline window with notebook and fix analysis window off-by-one"
```

---

### Task 8: Score spots during search aggregation

**Files:**
- Modify: `Methods/loadSlicesDMD.m:634-830` (search-aggregation block), `:556-561` and `:837-842` (Response map assembly), `:565-575` and `:845-855` (Stats assembly)

**Interfaces:**
- Consumes: `dmdFilterSweepsByRs`, `dmdHotspotVerdict`, `dmdSweepMeanStats`.
- Produces, on `responseMap`: `chargeExcessZ`, `chargeExcessPc`, `derivativeExcessZ`, `reanalysisHotspot`, `nSweepsKept` — each shaped like the existing `hotspot` field (per search, per depth, per spot). On `cellStats`: `filtered.max`, `filtered.min`, `filtered.maxTime`, `filtered.minTime`, `filtered.auc`, same shape.

- [ ] **Step 1: Drop Rs-failing sweeps before each spot is averaged**

After the depth loop fills `depthCurrentMap`, build a per-row keep mask from the depth's QC rows
(`Sweep`, `Rs`, `Rs_headerString`; `spotsAtDepth` carries `Sweep` per row, so the join is direct)
via `dmdFilterSweepsByRs(..., options.maxRs)`. Record `nSweepsKept` per spot. A spot with no
surviving sweep keeps an empty matrix rather than being removed.

- [ ] **Step 2: Build the reference per depth, falling back to the search**

Call `dmdHotspotVerdict('reference', ...)` on that depth's filtered spot traces. When the depth has
fewer than `criteria.minSpotsForReference` spots, rebuild it from every spot of the search.

- [ ] **Step 3: Score every spot and store the metrics**

`dmdHotspotVerdict('score', ...)` and `dmdHotspotVerdict('reanalysis', ...)` per spot, accumulated
into the five new Response-map fields alongside `searchHotspotMap`.

- [ ] **Step 4: Compute `cellStats.filtered.*`**

`dmdSweepMeanStats` on the same filtered traces over `[0, options.analysisWindowLength]`. Leave the
existing `cellStats.max/min/maxTime/minTime/auc` untouched.

- [ ] **Step 5: Verify on one real session**

Run `Shun_loadSliceData_DMD.m` sections 1–2 on one experiment folder. Confirm
`cells.("Response map"){1}.chargeExcessZ` is populated, finite for most spots, and that
`nSweepsKept` is at or below the repetition count for every spot.

- [ ] **Step 6: Commit**

```bash
git add Methods/loadSlicesDMD.m
git commit -m "feat: score DMD spots with notebook detection metrics at load time"
```

---

### Task 9: Consume the verdict in the search figures

**Files:**
- Modify: `Methods/analyzeDMDSearch.m:10-20` (options), `:252` (control slicing), `:272-276` (hotspot split), `:279-284` (stat vectors), `:347` (per-spot opacity)

**Interfaces:**
- Consumes: the Response-map and `cellStats.filtered` fields from Task 8.
- Produces: `options.hotspotType string = "evoked"` (`"acquisition"`/`"reanalysis"`/`"evoked"`), `options.minChargeZ double = 0.5`, `options.nullAlpha double = 0.3`.

- [ ] **Step 1: Replace the hotspot split with the verdict**

At `:273-274`, derive `hotspot_spotIdx` from `chargeExcessZ >= options.minChargeZ` for `"evoked"`,
`reanalysisHotspot` for `"reanalysis"`, and for `"acquisition"` from `isResponseMap` by looking up
each spot's box through `search_spotLocation`, as `loadSlicesDMD.m:2163` does. A `NaN` z is not a
hotspot. Error with a message naming the re-run if the new fields are absent.

- [ ] **Step 2: Point the statistics panels at the filtered stats**

At `:279-284`, read `curCell.Stats{1}.filtered.*` instead of `curCell.Stats{1}.*`.

- [ ] **Step 3: Drive per-spot opacity from the verdict**

At `:347`, `spotHotspot` comes from the same per-spot verdict rather than
`sum(depthHotspotFull{t})>=1`. The existing `opacity=1` / `opacity=options.nullAlpha` branch at
`:350-356` is otherwise unchanged. A spot with no surviving sweeps draws transparent. *(Review
Focus 1)*

- [ ] **Step 4: Fix the control-data slicing**

At `:252`, replace the fixed-length slice with the robust version commented at `:255-257`, which the
291-sample control window requires.

- [ ] **Step 5: Verify the figures on one real session**

Run `Shun_loadSliceData_DMD.m` section 3. Confirm a depth-1 search renders with no error, opaque
spots are those with `chargeExcessZ >= 0.5`, and `hotspotType="acquisition"` reproduces today's
`isResponseMap` spot for spot.

- [ ] **Step 6: Commit**

```bash
git add Methods/analyzeDMDSearch.m
git commit -m "feat: draw notebook hotspots opaque in DMD search figures"
```

---

### Task 10: Split the script's storage and plot windows

**Files:**
- Modify: `Shun_loadSliceData_DMD.m:18` (storage range), `:21-22`, `:64` (plot range)

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Set the storage window to `[-30,100]`**

In the `inputSessionParams_singleSlice` call at `:18`, `timeRange='[-30,100]'`. This is the value
passed to `loadSlices` and `loadSlicesDMD`.

- [ ] **Step 2: Pass the plot window to `analyzeSlice_DMD`**

At `:64`, pass `timeRange=[-10,50]` rather than the storage `timeRange`. `analyzeSlice_DMD`'s own
default is already `[-10,50]`, so no change is needed in that file.

- [ ] **Step 3: Verify end to end**

Re-run the script from a clean `Results-*` folder. Confirm stored traces span −30 to +100 ms, the
figure tiles still span −10 to +50 ms, and the notebook's `charge_excess_z` matches MATLAB's
`chargeExcessZ` for one cell, search and depth.

- [ ] **Step 4: Commit**

```bash
git add Shun_loadSliceData_DMD.m
git commit -m "feat: store a wider DMD window than the figures display"
```
