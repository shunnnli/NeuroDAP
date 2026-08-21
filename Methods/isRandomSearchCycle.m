function tf = isRandomSearchCycle(cycleName)

% Return true for random-search cycle names independent of formatting.
% Examples matched: randomSearch, Random_Search, and random search.cyc.

normalizedName = regexprep(lower(string(cycleName)), '[^a-z0-9]', '');
tf = contains(normalizedName, "randomsearch");

end
