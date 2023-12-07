function stats = analyzeStages(trace,stageTime,options)

arguments
    trace double
    stageTime double
    
    options.timeRange double = [-15,15]
    options.nboot double = 10000
    options.finalFs double = 50
    options.lick_binSize double
end

% Provide lick_binSize if analyzing lick traces
if isfield(options,'lick_binSize'); options.finalFs = 1/options.lick_binSize; end

% Calculate subtrial averages
stageAvg = nan(size(trace,1),size(stageTime,1));
stageBin = (stageTime - options.timeRange(1)) * options.finalFs;
stageAvgFit = nan(size(stageTime,1),2);
stageAvgFit_bs = nan(size(stageTime,1),options.nboot,2);
pAvg_intercept = nan(size(stageTime,1),2);
pAvg_slope = nan(size(stageTime,1),2);
for stage = 1:size(stageTime,1)
    stageWindow = stageBin(stage,1):stageBin(stage,2);
    stageAvg(:,stage) = mean(trace(:,stageWindow),2);

    % Fit stageAvg across session
    x = 1:size(trace,1); y = stageAvg(:,stage)';
    if length(y) <= 1; continue; end
    stageAvgFit(stage,:) = polyfit(x,y,1);

    % Test for significance (shuffled individual points, 
    % calculate slope to build distribution, and test significance)
    for sample = 1:options.nboot
        bs_data = y(randperm(length(y)));
        stageAvgFit_bs(stage,sample,:) = polyfit(x,bs_data,1);
    end
    % Calculate one-side p value
    pAvg_intercept(stage,1) = sum(stageAvgFit_bs(stage,:,end)<=stageAvgFit(stage,end))/options.nboot;
    pAvg_intercept(stage,2) = sum(stageAvgFit_bs(stage,:,end)>=stageAvgFit(stage,end))/options.nboot;
    pAvg_slope(stage,1) = sum(stageAvgFit_bs(stage,:,end-1)<=stageAvgFit(stage,end-1))/options.nboot;
    pAvg_slope(stage,2) = sum(stageAvgFit_bs(stage,:,end-1)>=stageAvgFit(stage,end-1))/options.nboot;
end

% Calculate subtrial peaks
stageMax = nan(size(trace,1),size(stageTime,1));
stageMaxFit = nan(size(stageTime,1),2);
stageMaxFit_bs = nan(size(stageTime,1),options.nboot,2);
pMax_intercept = nan(size(stageTime,1),2);
pMax_slope = nan(size(stageTime,1),2);
for stage = 1:size(stageTime,1)
    stageWindow = stageBin(stage,1):stageBin(stage,2);
    stageMax(:,stage) = max(trace(:,stageWindow),[],2);

    % Fit stageMax across session
    x = 1:size(trace,1); y = stageMax(:,stage)';
    if length(y) <= 1; continue; end
    stageMaxFit(stage,:) = polyfit(x,y,1);

    % Test for significance (shuffled individual points, 
    % calculate slope to build distribution, and test significance)
    for sample = 1:options.nboot
        bs_data = y(randperm(length(y)));
        stageMaxFit_bs(stage,sample,:) = polyfit(x,bs_data,1);
    end
    % Calculate one-side p value
    pMax_intercept(stage,1) = sum(stageMaxFit_bs(stage,:,end)<=stageMaxFit(stage,end))/options.nboot;
    pMax_intercept(stage,2) = sum(stageMaxFit_bs(stage,:,end)>=stageMaxFit(stage,end))/options.nboot;
    pMax_slope(stage,1) = sum(stageMaxFit_bs(stage,:,end-1)<=stageMaxFit(stage,end-1))/options.nboot;
    pMax_slope(stage,2) = sum(stageMaxFit_bs(stage,:,end-1)>=stageMaxFit(stage,end-1))/options.nboot;
end

% Calculate subtrial troughs
stageMin = nan(size(trace,1),size(stageTime,1));
stageMinFit = nan(size(stageTime,1),2);
stageMinFit_bs = nan(size(stageTime,1),options.nboot,2);
pMin_intercept = nan(size(stageTime,1),2);
pMin_slope = nan(size(stageTime,1),2);
for stage = 1:size(stageTime,1)
    stageWindow = stageBin(stage,1):stageBin(stage,2);
    stageMin(:,stage) = min(trace(:,stageWindow),[],2);

    % Fit stageMin across session
    x = 1:size(trace,1); y = stageMin(:,stage)';
    if length(y) <= 1; continue; end
    stageMinFit(stage,:) = polyfit(x,y,1);

    % Test for significance (shuffled individual points, 
    % calculate slope to build distribution, and test significance)
    for sample = 1:options.nboot
        bs_data = y(randperm(length(y)));
        stageMinFit_bs(stage,sample,:) = polyfit(x,bs_data,1);
    end
    % Calculate one-side p value
    pMin_intercept(stage,1) = sum(stageMinFit_bs(stage,:,end)<=stageMinFit(stage,end))/options.nboot;
    pMin_intercept(stage,2) = sum(stageMinFit_bs(stage,:,end)>=stageMinFit(stage,end))/options.nboot;
    pMin_slope(stage,1) = sum(stageMinFit_bs(stage,:,end-1)<=stageMinFit(stage,end-1))/options.nboot;
    pMin_slope(stage,2) = sum(stageMinFit_bs(stage,:,end-1)>=stageMinFit(stage,end-1))/options.nboot;
end

% Store in stats
stats.stageAvg.data = stageAvg;
stats.stageAvg.fit = stageAvgFit;
stats.stageAvg.stats.bs = stageAvgFit_bs;
stats.stageAvg.stats.pval_intercept = pAvg_intercept;
stats.stageAvg.stats.pval_slope = pAvg_slope;

stats.stageMax.data = stageMax;
stats.stageMax.fit = stageMaxFit;
stats.stageMax.stats.bs = stageMaxFit_bs;
stats.stageMax.stats.pval_intercept = pMax_intercept;
stats.stageMax.stats.pval_slope = pMax_slope;

stats.stageMin.data = stageMin;
stats.stageMin.fit = stageMinFit;
stats.stageMin.stats.bs = stageMinFit_bs;
stats.stageMin.stats.pval_intercept = pMin_intercept;
stats.stageMin.stats.pval_slope = pMin_slope;

stats.stageAvg.stageTime = stageTime;
stats.stageMax.stageTime = stageTime;
stats.stageMin.stageTime = stageTime;

stats.options = options;

end