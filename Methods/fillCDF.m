function filled = fillCDF(cdf)
% fillCDF  Remove empty rows and fill NaNs in each row of an interpolated CDF matrix
%
%   filled = fillCDF(cdf)
%
%   Input:
%     cdf       – m×n matrix, each row is an ECDF on a common grid, possibly containing NaNs
%   Output:
%     filled    – p×n matrix (p ≤ m), with:
%                   • rows that were all-NaN removed
%                   • leading NaNs → 0
%                   • trailing NaNs → 1
%                   • intermediate NaNs → carried forward

    % 1) drop rows that are all NaN
    keep = ~all(isnan(cdf), 2);
    cdf  = cdf(keep, :);
    
    [m, n] = size(cdf);
    filled = zeros(m, n);

    for i = 1:m
        row = cdf(i, :);
        idx = find(~isnan(row));

        if isempty(idx)
            % (shouldn’t happen now, since we removed all-NaN rows)
            row_filled = zeros(1, n);
        else
            firstValid = idx(1);
            lastValid  = idx(end);

            % 1) before firstValid → 0
            row(1:firstValid-1) = 0;
            % 2) after  lastValid  → 1
            row(lastValid+1:end) = 1;
            % 3) carry-forward through any remaining NaNs
            row_filled = fillmissing(row, 'previous', 2);
        end

        filled(i, :) = row_filled;
    end
end

