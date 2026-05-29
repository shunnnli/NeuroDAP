function amp = nearestProtocolAmp(value, expected_amps, tol)
    [d, idx] = min(abs(expected_amps - value));
    if d <= tol
        amp = expected_amps(idx);
    else
        amp = round(value, 3);
    end
end