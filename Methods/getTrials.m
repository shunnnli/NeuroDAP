function [allTrials,pairIdx] = getTrials(varargin,options)

arguments (Repeating)
    varargin
end

arguments
    options.threshold double = 2.5
    options.behaviorFs double = 10000;
end

% Add stim and tone together
combined = sort(horzcat(varargin{:}));
diff_combined = [5*options.behaviorFs, diff(combined)];

if isempty(combined); allTrials = []; return; end

% Combine stim and tone that are sufficiently close to each other
allTrials = combined(diff_combined > options.threshold*options.behaviorFs);
pairIdx = combined(diff_combined < options.threshold*options.behaviorFs);

end