function [amplitudes,peakIndices,minIndices] = detectPaperPositivePeaks(x,minimumRise)
% Positive local-minimum-to-next-maximum rises in percent dF/F.
% McGregor et al., doi:10.64898/2026.07.24.740623, Fiber photometry.
% The paper specifies >=1.5% rise, not prominence, absolute height or a
% minimum distance. Pairing with the immediately preceding local minimum
% and taking the first sample of a flat extremum are implementation choices.
% Endpoints are not extrema; never bridge a missing sample or trial boundary.
x = double(x(:)');
assert(all(isfinite(x)),'Pass separate finite segments to the detector.');
assert(isscalar(minimumRise) && minimumRise>0,'Minimum rise must be positive.');
% Collapse plateaus to avoid detecting multiple peaks on a flat top.
idx = [1,find(diff(x)~=0)+1]; y = x(idx);
turn = diff(sign(diff(y)));
mins = find(turn>0)+1; maxs = find(turn<0)+1;
amplitudes = []; peakIndices = []; minIndices = [];
for k = maxs
    m = mins(find(mins<k,1,'last'));
    if isempty(m); continue; end
    rise = y(k)-y(m);
    if rise>=minimumRise
        amplitudes(end+1) = rise; %#ok<AGROW>
        peakIndices(end+1) = idx(k); %#ok<AGROW>
        minIndices(end+1) = idx(m); %#ok<AGROW>
    end
end
end
