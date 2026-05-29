function events = findRampEvents(target_dff, fs, expected_amps, threshold, amp_tol)
    slope_threshold_per_s = 0.015;

    y = target_dff(:);
    smooth_win = max(3, round(0.05 * fs));
    y_smooth = movmedian(y, smooth_win, "omitnan");

    slope = gradient(y_smooth) * fs;
    slope_mask = abs(slope) > slope_threshold_per_s;

    segments = contiguousRegions(slope_mask);
    segments = mergeCloseSegments(segments, round(0.30 * fs));

    min_samples = round(0.5 * fs);
    max_samples = round(4.0 * fs);
    level_win = max(1, round(0.25 * fs));

    rows = [];

    for i = 1:size(segments, 1)
        start_idx = segments(i, 1);
        end_idx = segments(i, 2);
        n_samples = end_idx - start_idx + 1;

        if n_samples < min_samples || n_samples > max_samples
            continue
        end

        pre_idx = max(1, start_idx - level_win):(start_idx - 1);
        post_idx = (end_idx + 1):min(numel(y), end_idx + level_win);

        if isempty(pre_idx)
            start_level = y_smooth(start_idx);
        else
            start_level = median(y(pre_idx), "omitnan");
        end

        if isempty(post_idx)
            end_level = y_smooth(end_idx);
        else
            end_level = median(y(post_idx), "omitnan");
        end

        delta = end_level - start_level;

        if abs(delta) < threshold
            continue
        end

        if abs(start_level) <= abs(end_level)
            ramp_phase = "away_from_baseline";
            amp_raw = end_level;
        else
            ramp_phase = "return_to_baseline";
            amp_raw = start_level;
        end

        amp = nearestProtocolAmp(amp_raw, expected_amps, amp_tol);

        rows = [rows; struct( ...
            "session_type", "ramp", ...
            "kind", "ramp", ...
            "ramp_phase", ramp_phase, ...
            "trial", floor(numel(rows) / 2) + 1, ...
            "amp", amp, ...
            "amp_raw", amp_raw, ...
            "start_level", start_level, ...
            "end_level", end_level, ...
            "delta", delta, ...
            "start_idx", start_idx, ...
            "end_idx", end_idx, ...
            "start_s", (start_idx - 1) / fs, ...
            "end_s", (end_idx - 1) / fs, ...
            "duration_s", n_samples / fs)];
    end

    events = struct2table(rows);
end


function regions = contiguousRegions(mask)
    mask = logical(mask(:));
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    regions = [starts, ends];
end

function merged = mergeCloseSegments(segments, max_gap_samples)
    if isempty(segments)
        merged = segments;
        return
    end

    merged = segments(1, :);

    for i = 2:size(segments, 1)
        start_idx = segments(i, 1);
        end_idx = segments(i, 2);

        if start_idx - merged(end, 2) <= max_gap_samples
            merged(end, 2) = end_idx;
        else
            merged = [merged; start_idx, end_idx];
        end
    end
end