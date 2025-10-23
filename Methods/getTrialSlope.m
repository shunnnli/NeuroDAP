function slope_amp = getTrialSlope(trial_amp, trial_window, options)
% slope_amp(:,t): LS slope over trials t+trial_window(1) : t+trial_window(2), per row.
% options.edge:
%   'omit' -> require the FULL window in-bounds with no NaNs; otherwise NaN
%   'pad'  -> use ONLY the in-bounds part of the window (>=2 pts), then
%             time-wise interpolate remaining gaps
% NOTE: if the current trial A(:,t) is NaN, slope_amp(:,t) = NaN (always).

arguments
    trial_amp double
    trial_window (1,2) double = [-5 5]
    options.edge (1,1) string {mustBeMember(options.edge,["omit","pad"])} = "omit"
end

A = trial_amp; wasRow = isvector(A) && isrow(A); wasCol = isvector(A) && iscolumn(A);
if isvector(A), A = A(:)'; end
[M,T] = size(A);
slope_amp = nan(M,T);
curOK = ~isnan(A);                         % must be valid at current trial

% integer window offsets and x-coordinates (in "trials")
offsets = round(trial_window(1)) : round(trial_window(2));
x_full  = offsets(:);                      % column
L = numel(x_full);

switch lower(options.edge)
    case 'omit'    % need the entire window in-bounds
        Sx  = sum(x_full); 
        Sxx = sum(x_full.^2); 
        denom = L*Sxx - Sx^2;
        for t = 1:T
            idx = t + offsets;
            if all(idx>=1 & idx<=T)
                W   = A(:,idx);                 % M x L
                Sy  = sum(W,2); 
                Sxy = W * x_full;              % (M x 1)
                s = (L*Sxy - Sx*Sy) / denom;
                s(~curOK(:,t)) = NaN;          % force NaN if current trial is NaN
                slope_amp(:,t) = s;
            end
        end

    otherwise      % 'pad' -> use in-bounds subset (>=2 pts)
        for t = 1:T
            idx = t + offsets;
            keep = idx>=1 & idx<=T;
            k = nnz(keep);
            if k >= 2
                x   = x_full(keep);
                W   = A(:, idx(keep));         % M x k
                Sx  = sum(x); 
                Sxx = sum(x.^2); 
                denom = k*Sxx - Sx^2;
                Sy  = sum(W,2); 
                Sxy = W * x;                   % (M x 1)
                s = (k*Sxy - Sx*Sy) / denom;   % rows with NaNs in W -> NaN
                s(~curOK(:,t)) = NaN;          % do not compute if current trial is NaN
                slope_amp(:,t) = s;
            end
        end
        % smooth over time-wise gaps, but never fill current-NaN trials
        slope_amp = fillmissing(slope_amp,'linear',2);
        slope_amp = fillmissing(slope_amp,'next',  2);
end

% enforce current-NaN rule after any filling
slope_amp(~curOK) = NaN;

% restore original vector orientation
if wasCol, slope_amp = slope_amp.'; end
if wasRow, slope_amp = slope_amp;   end
end
