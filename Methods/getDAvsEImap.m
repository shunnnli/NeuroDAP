function DAvsEImap = getDAvsEImap(DAtrend,animalEIindex,options)

% Get DA vs EI map
arguments
    DAtrend struct
    animalEIindex double

    options.metric string = 'cosine'

    options.mapType string = 'slope'
    options.nTrials double = 30

    % P value options
    options.pval logical = true
    options.nboot double = 1000
end

%% Check input

if sum(contains(options.metric,["cosine","cos"]))
    metric = 'cosine';
elseif contains(options.metric,"slope")
    metric = 'slope';
end


if sum(contains(options.mapType,["slope","slopes"]))
    mapType = 'slopeMap';
elseif sum(contains(options.mapType,["diff","dif"]))
    mapType = 'diffMap';
end

%% Setup

fields = {'max', 'min', 'avg', 'amp'};
nAnimals = length(unique({DAtrend.animal}));

% Preallocate map
DAvsEImap = struct();
for f = 1:length(fields)
    DAvsEImap.(fields{f}).raw = nan(options.nTrials*2,options.nTrials*2);
    DAvsEImap.(fields{f}).smoothed = nan(options.nTrials*2,options.nTrials*2);

    if options.pval
        DAvsEImap.(fields{f}).pval_raw = nan(options.nTrials*2,options.nTrials*2);
        DAvsEImap.(fields{f}).pval_smoothed = nan(options.nTrials*2,options.nTrials*2);
    end
end

%% Generate trialIdx

trialIdx = zeros(nAnimals,options.nTrials*2);
for a = 1:nAnimals
    maxTrials = DAtrend(a).nTrials;
    forwardIdx = 1:options.nTrials;
    backwardIdx = flip(maxTrials - forwardIdx + 1);
    trialIdx(a,:) = [forwardIdx,backwardIdx];
end

%% Generate DA vs EI map

disp('Ongoing: generating DA vs EI map');
for i = 1:size(trialIdx,2)
    for j = i+1:size(trialIdx,2)

        % Get raw DA stat
        DAstat_raw.max = nan(nAnimals,1);
        DAstat_raw.min = nan(nAnimals,1);
        DAstat_raw.avg = nan(nAnimals,1);
        DAstat_raw.amp = nan(nAnimals,1);
        for a = 1:nAnimals
            t1 = trialIdx(a,i);
            t2 = trialIdx(a,j);
            animalDAStat = getDAtrend(DAtrend,t1,t2,animalIdx=a,mapType=mapType,dataType='raw');
            DAstat_raw.max(a) = animalDAStat.max;
            DAstat_raw.min(a) = animalDAStat.min;
            DAstat_raw.avg(a) = animalDAStat.avg;
            DAstat_raw.amp(a) = animalDAStat.amp;
        end

        % Get smooth DA stat
        DAstat_smoothed.max = nan(nAnimals,1);
        DAstat_smoothed.min = nan(nAnimals,1);
        DAstat_smoothed.avg = nan(nAnimals,1);
        DAstat_smoothed.amp = nan(nAnimals,1);
        for a = 1:nAnimals
            t1 = trialIdx(a,i);
            t2 = trialIdx(a,j);
            animalDAStat = getDAtrend(DAtrend,t1,t2,animalIdx=a,mapType=mapType,dataType='smooth');
            DAstat_smoothed.max(a) = animalDAStat.max;
            DAstat_smoothed.min(a) = animalDAStat.min;
            DAstat_smoothed.avg(a) = animalDAStat.avg;
            DAstat_smoothed.amp(a) = animalDAStat.amp;
        end

        % Skip if animal number is less than 5
        validIdx = ~isnan(DAstat_raw.max);
        if sum(validIdx) <= 5; continue; end

        % Calculate cosine similarity for DA vs EI
        for f = 1:length(fields)

            % Get raw heatmap
            X = DAstat_raw.(fields{f})(validIdx); 
            Y = animalEIindex(validIdx);
            if strcmpi(options.metric,'cosine')
                metric_obs = dot(X, Y) / (norm(X) * norm(Y));
                % Calculate p value
                cosSimFun = @(idx) dot(X(idx), Y(idx)) / (norm(X(idx)) * norm(Y(idx)));
                boot_cosSim = bootstrp(options.nboot, cosSimFun, 1:length(Y));
                metric_pval = sum(abs(boot_cosSim) >= abs(metric_obs))/options.nboot;
            elseif strcmpi(options.metric,'slope')
                fit_obs = polyfit(X, Y, 1); metric_obs = fit_obs(1);
                % Calculate p value
                fit_boot = bootstrp(options.nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
                metric_pval = sum(abs(fit_boot(:,1)) >= abs(metric_obs))/options.nboot;
            end
            DAvsEImap.(fields{f}).raw(i,j) = metric_obs;
            DAvsEImap.(fields{f}).pval_raw(i,j) = metric_pval;

            % Get smoothed heatmap
            X = DAstat_smoothed.(fields{f})(validIdx); 
            Y = animalEIindex(validIdx);
            if strcmpi(options.metric,'cosine')
                metric_obs = dot(X, Y) / (norm(X) * norm(Y));
                % Calculate p value
                cosSimFun = @(idx) dot(X(idx), Y(idx)) / (norm(X(idx)) * norm(Y(idx)));
                boot_cosSim = bootstrp(options.nboot, cosSimFun, 1:length(Y));
                metric_pval = sum(abs(boot_cosSim) >= abs(metric_obs))/options.nboot;
            elseif strcmpi(options.metric,'slope')
                fit_obs = polyfit(X, Y, 1); metric_obs = fit_obs(1);
                % Calculate p value
                fit_boot = bootstrp(options.nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
                metric_pval = sum(abs(fit_boot(:,1)) >= abs(metric_obs))/options.nboot;
            end
            DAvsEImap.(fields{f}).smoothed(i,j) = metric_obs;
            DAvsEImap.(fields{f}).pval_smoothed(i,j) = metric_pval;
        end
    end
end

DAvsEImap.options = options;
disp('Finished: DA vs EI map generated');

end