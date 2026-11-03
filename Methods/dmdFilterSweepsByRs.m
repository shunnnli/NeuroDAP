function keep = dmdFilterSweepsByRs(rs,rsHeader,maxRs)

% Which sweeps survive the series-resistance limit.
%
% The RC-step fit fails outright on a minority of sweeps, where Rs reaches
% several hundred MOhm while the header reads tens, and the header is written
% as a flat 0 on others. Judging the minimum of the two reads each sweep as
% "the series resistance was at most this bad": it keeps sweeps a failed fit
% would have thrown away, and drops sweeps the header alone would have kept.
%
% A non-positive header value is not a measurement, so it is treated as
% missing rather than as an unbeatably small Rs. A sweep with no usable
% estimate at all does not pass.
%
% Mirrors resolve_qc_metrics(choice="min") in pyNeuroDAP's Shun_dmd_analysis
% notebook, whose SWEEP_QC_LIMITS is {"Rs": 30.0}.

arguments
    rs double
    rsHeader double
    maxRs (1,1) double
end

rsHeader(rsHeader <= 0) = NaN;

resolved = min([rs(:), rsHeader(:)],[],2,'omitnan');
keep = reshape(resolved <= maxRs, size(rs));

end
