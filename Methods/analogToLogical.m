function onVec = analogToLogical(x, Fs, options)

arguments
    x double
    Fs double = 10000

    options.noiseThresholdPct double = 1   % threshold as % of signal range above OFF level
    options.smooth_ms double = 3
end

    x = x(:);

    % 1) Smooth slightly to suppress small analog noise
    win = max(1, round(options.smooth_ms/1000 * Fs));
    xSmooth = movmedian(x, win);

    % 2) Estimate OFF and ON levels robustly from the signal itself
    xOff = prctile(xSmooth, 5);
    xOn  = prctile(xSmooth, 95);

    % If signal has almost no dynamic range, return all false
    if (xOn - xOff) < eps
        onVec = false(size(x));
        return
    end

    % 3) Single threshold: samples below threshold are OFF, others are ON
    thr = xOff + (options.noiseThresholdPct/100) * (xOn - xOff);
    onVec = xSmooth >= thr;

    % Match input shape
    onVec = reshape(onVec, size(x));
end