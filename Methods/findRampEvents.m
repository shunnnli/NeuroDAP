function events = findRampEvents(target_dff, fs, expected_amps, threshold, amp_tol)
    active = abs(target_dff) > threshold;
    segments = contiguousRegions(active);
    segments = mergeCloseSegments(segments, round(0.15 * fs));

    min_samples = round(1.0 * fs);
    trim = round(0.2 * fs);
    rows = [];

    for i = 1:size(segments, 1)
        start_idx = segments(i, 1);
        end_idx = segments(i, 2);

        if end_idx - start_idx + 1 < min_samples
            continue
        end

        lo = min(end_idx, start_idx + trim);
        hi = max(lo, end_idx - trim);
        inner = target_dff(lo:hi);

        amp_raw = median(inner, "omitnan");
        amp = nearestProtocolAmp(amp_raw, expected_amps, amp_tol);

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