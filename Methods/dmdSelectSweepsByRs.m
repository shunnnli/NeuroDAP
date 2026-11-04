function keep = dmdSelectSweepsByRs(rs,rsHeader,maxRs,fallbackMaxRs)
% Relax the Rs limit only for a spot with exactly one initially usable sweep.
arguments
    rs double
    rsHeader double
    maxRs (1,1) double
    fallbackMaxRs (1,1) double = 35
end
keep = dmdFilterSweepsByRs(rs,rsHeader,maxRs);
if nnz(keep) == 1 && fallbackMaxRs > maxRs
    keep = dmdFilterSweepsByRs(rs,rsHeader,fallbackMaxRs);
end
end
