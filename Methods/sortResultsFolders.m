function resultsList = sortResultsFolders(parentPath)

% Results-* folders under parentPath, OLDEST FIRST, so {end,:} is the newest.
%
% Returns the same N x 6 cell array struct2cell(dir(...))' gives -- name,
% folder, date, bytes, isdir, datenum -- so callers index it as before.
%
% Sorted on the folder NAME, not on dir's date. Two reasons:
%
%   * dir's date is a char vector like '28-Oct-2026 17:43:21'. Sorted as text,
%     '03-Nov-2026' lands BEFORE '28-Oct-2026', so {end} returned the older
%     folder whenever the newer one fell on a smaller day of the month. A
%     freshly written Results folder was then silently skipped in favour of a
%     stale one.
%   * a folder's modification time changes whenever anything inside it is
%     written, so even a correct chronological sort on it would not mean
%     "the results from the latest run".
%
% The names are Results-YYYYMMDD, zero padded, so sorting them as text is
% chronological. It also matches the sibling branch that selects a folder by
% date, which already compares against the name.

arguments
    parentPath string
end

listing = dir(fullfile(parentPath,'Results-*'));
listing = listing([listing.isdir]);

if isempty(listing)
    resultsList = {};
    return
end

resultsList = sortrows(struct2cell(listing)',1);

end
