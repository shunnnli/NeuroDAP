function events = findStepEvents(target_dff, fs, expected_amps, zero_threshold, amp_tol)

    y = target_dff(:);

    % Smooth a little so noise does not create fake edges.
    smooth_win = max(3, round(0.05 * fs));
    y_smooth = movmedian(y, smooth_win, "omitnan");

    % Detect transitions between command levels.
    step_change_threshold = zero_threshold;
    edge_idx = find(abs(diff(y_smooth)) > step_change_threshold) + 1;

    % Merge edge samples that belong to the same transition.
    if ~isempty(edge_idx)
        edge_regions = contiguousRegionsFromIdx(edge_idx, round(0.15 * fs));
        change_idx = round(mean(edge_regions, 2));
    else
        change_idx = [];
    end

    % Plateaus are the regions between transitions.
    bounds = [1; change_idx(:); numel(y)];
    rows = [];
    min_samples = round(0.5 * fs);
    edge_pad = round(0.15 * fs);  % ignore transition artifact at plateau edges

    for i = 1:(numel(bounds) - 1)
        start_idx = bounds(i);
        end_idx = bounds(i + 1) - 1;

        if end_idx - start_idx + 1 < min_samples
            continue
        end

        lo = min(end_idx, start_idx + edge_pad);
        hi = max(lo, end_idx - edge_pad);

        amp_raw = median(y(lo:hi), "omitnan");

        % Only keep plateaus that are meaningfully away from 0.
        if abs(amp_raw) < zero_threshold
            continue
        end

        amp = nearestProtocolAmp(amp_raw, expected_amps, amp_tol);

        % Keep amp tolerance: skip plateaus too far from provided amps.
        if isnan(amp)
            continue
        end

        rows = [rows; struct(...
            "kind", "step", ...
            "trial", numel(rows) + 1, ...
            "amp", amp, ...
            "amp_raw", amp_raw, ...
            "start_idx", start_idx, ...
            "end_idx", end_idx, ...
            "start_s", (start_idx - 1) / fs, ...
            "end_s", (end_idx - 1) / fs, ...
            "duration_s", (end_idx - start_idx + 1) / fs)];
    end

    if isempty(rows)
        events = table();
    else
        events = struct2table(rows);
    end
end


function regions = contiguousRegionsFromIdx(idx, max_gap_samples)
    idx = idx(:);

    if isempty(idx)
        regions = [];
        return
    end

    starts = idx([true; diff(idx) > max_gap_samples]);
    ends = idx([diff(idx) > max_gap_samples; true]);

    regions = [starts, ends];
end