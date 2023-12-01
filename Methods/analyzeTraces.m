function analysis = analyzeTraces(timeSeries,lick,analysisEvents,analysisLabels,params,options)

arguments
    timeSeries struct
    lick double
    analysisEvents cell
    analysisLabels cell
    params

    options.task string = "NaN"
    options.timeRange double = [-15, 15]
    options.stageTime double = [-2,0;0,2]
    options.lick_binSize double = 0.1
    options.nboot double = 10000
    
    options.trialTable table
    options.trialNumber cell

end

disp('Ongoing: analyze traces and saved in anlaysis struct');
options.analysisDate = char(datetime('today','Format','yyyyMMdd'));

%% Define analysis params
analysis = struct([]);
timeRange = options.timeRange;
stageTime = options.stageTime;

% Define tasks and other properties
if isfield(params.session,'task')
    options.task = params.session.task;
end

% Make sure trialTable and trialNumber is provided
if ~isfield(options,'trialNumber') || ~isfield(options,'trialTable')
    warning('analyzeTraces: Strongly recommend to include trialTable and trialNumber!')
    disp('Ongoing: starting a 30s timer, stop and modify the code if needed!');
    pause(30);
    disp('Finished: 30s timer passed, code resume as is');
end

% Change all trialNumber to nx1 vector
for i = 1:length(options.trialNumber)
    options.trialNumber{i} = reshape(options.trialNumber{i},[length(options.trialNumber{i}),1]);
end

%% Loop through all events
for i = 1:length(analysisEvents)
    if isempty(analysisEvents{i})
        disp(['     Ongoing: ',analysisLabels{i},' is empty, skipped!']);
        continue
    end
    %% 1. Loop through all timeSeries
    for signal = 1:size(timeSeries,2)
        disp(['     Ongoing: analyze ',analysisLabels{i},' in ',timeSeries(signal).name,' signal']);

        % Load signal of interest
        data = timeSeries(signal).data;
        finalFs = round(timeSeries(signal).finalFs);
        system = timeSeries(signal).system;

        % Save overall traces
        [trace,t] = plotTraces(analysisEvents{i},timeRange,data,[1,1,1],params,...
                        signalFs=finalFs,signalSystem=system,plot=false);

        % Calculate subtrial averages
        stageAvg = nan(size(trace,1),size(stageTime,1));
        stageBin = (stageTime - timeRange(1)) * finalFs;
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

        % Save traces and anlaysis data
        row = size(analysis,2) + 1;
        analysis(row).animal = params.session.animal;
        analysis(row).date = params.session.date;
        analysis(row).session = params.session.name;
        analysis(row).task = options.task;
        analysis(row).event = analysisLabels{i};
        analysis(row).name = timeSeries(signal).name;
        analysis(row).system = system;
        analysis(row).data = trace;
        analysis(row).timestamp = t;
        analysis(row).timeRange = timeRange;
        analysis(row).finalFs = finalFs;

        analysis(row).stageAvg.data = stageAvg;
        analysis(row).stageAvg.fit = stageAvgFit;
        analysis(row).stageAvg.stats.bs = stageAvgFit_bs;
        analysis(row).stageAvg.stats.pval_intercept = pAvg_intercept;
        analysis(row).stageAvg.stats.pval_slope = pAvg_slope;

        analysis(row).stageMax.data = stageMax;
        analysis(row).stageMax.fit = stageMaxFit;
        analysis(row).stageMax.stats.bs = stageMaxFit_bs;
        analysis(row).stageMax.stats.pval_intercept = pMax_intercept;
        analysis(row).stageMax.stats.pval_slope = pMax_slope;

        analysis(row).stageMin.data = stageMin;
        analysis(row).stageMin.fit = stageMinFit;
        analysis(row).stageMin.stats.bs = stageMinFit_bs;
        analysis(row).stageMin.stats.pval_intercept = pMin_intercept;
        analysis(row).stageMin.stats.pval_slope = pMin_slope;

        analysis(row).stageAvg.stageTime = stageTime;
        analysis(row).stageMax.stageTime = stageTime;
        analysis(row).stageMin.stageTime = stageTime;
        
        analysis(row).trialInfo.trialNumber = options.trialNumber{i};
        if ~isempty(options.trialTable)
            analysis(row).trialInfo.trialTable = options.trialTable(options.trialNumber{i},:);
        end
        analysis(row).options = options;
    end
    
    %% 2. Get lick events
    disp(['     Ongoing: analyze ',analysisLabels{i},' in lick signal']);
    [lickRate,lickEvents,t_licks] = plotLicks(analysisEvents{i},timeRange,...
                    options.lick_binSize,[1 1 1],[],lick,params,plot=false);

    % Calculate subtrial averages
    stageAvg_lick = nan(size(lickRate,1),size(stageTime,1));
    stageBin = (stageTime - timeRange(1)) * (1/options.lick_binSize);
    stageAvgFit_lick = nan(size(stageTime,1),2);
    stageAvgFit_bs = nan(size(stageTime,1),options.nboot,2);
    pAvg_intercept = nan(size(stageTime,1),2);
    pAvg_slope = nan(size(stageTime,1),2);
    for stage = 1:size(stageTime,1)
        stageWindow = stageBin(stage,1):stageBin(stage,2);
        stageAvg_lick(:,stage) = mean(lickRate(:,stageWindow),2);

        % Fit stageAvg across session
        x = 1:size(lickRate,1); y = stageAvg_lick(:,stage)';
        if length(y) <= 1; continue; end
        stageAvgFit_lick(stage,:) = polyfit(x,y,1);

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
    stageMax_lick = nan(size(lickRate,1),size(stageTime,1));
    stageMaxFit_lick = nan(size(stageTime,1),2);
    stageMaxFit_bs = nan(size(stageTime,1),options.nboot,2);
    pMax_intercept = nan(size(stageTime,1),2);
    pMax_slope = nan(size(stageTime,1),2);
    for stage = 1:size(stageTime,1)
        stageWindow = stageBin(stage,1):stageBin(stage,2);
        stageMax_lick(:,stage) = max(lickRate(:,stageWindow),[],2);

        % Fit stage peak across session
        x = 1:size(lickRate,1); y = stageMax_lick(:,stage)';
        if length(y) <= 1; continue; end
        stageMaxFit_lick(stage,:) = polyfit(x,y,1);

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
    stageMin_lick = nan(size(lickRate,1),size(stageTime,1));
    stageMinFit_lick = nan(size(stageTime,1),2);
    stageMinFit_bs = nan(size(stageTime,1),options.nboot,2);
    pMin_intercept = nan(size(stageTime,1),2);
    pMin_slope = nan(size(stageTime,1),2);
    for stage = 1:size(stageTime,1)
        stageWindow = stageBin(stage,1):stageBin(stage,2);
        stageMin_lick(:,stage) = min(lickRate(:,stageWindow),[],2);

        % Fit stage trough across session
        x = 1:size(lickRate,1); y = stageMin_lick(:,stage)';
        if length(y) <= 1; continue; end
        stageMinFit_lick(stage,:) = polyfit(x,y,1);

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

    % Save lick traces
    row = size(analysis,2) + 1;
    analysis(row).animal = params.session.animal;
    analysis(row).date = params.session.date;
    analysis(row).session = params.session.name;
    analysis(row).task = options.task;
    analysis(row).event = analysisLabels{i};
    analysis(row).name = 'Lick';
    analysis(row).system = 'Lick';
    analysis(row).data.lickRate = lickRate;
    analysis(row).data.lickEvents = lickEvents;
    analysis(row).timestamp = t_licks;
    analysis(row).timeRange = timeRange;
    analysis(row).finalFs = 1/options.lick_binSize;

    analysis(row).stageAvg.data = stageAvg_lick;
    analysis(row).stageAvg.fit = stageAvgFit_lick;
    analysis(row).stageAvg.stats.bs = stageMaxFit_bs;
    analysis(row).stageAvg.stats.pval_intercept = pAvg_intercept;
    analysis(row).stageAvg.stats.pval_slope = pAvg_slope;
    
    analysis(row).stageMax.data = stageMax_lick;
    analysis(row).stageMax.fit = stageMaxFit_lick;
    analysis(row).stageMax.stats.bs = stageMaxFit_bs;
    analysis(row).stageMax.stats.pval_intercept = pMax_intercept;
    analysis(row).stageMax.stats.pval_slope = pMax_slope;

    analysis(row).stageMin.data = stageMin_lick;
    analysis(row).stageMin.fit = stageMinFit_lick;
    analysis(row).stageMin.stats.bs = stageMinFit_bs;
    analysis(row).stageMin.stats.pval_intercept = pMin_intercept;
    analysis(row).stageMin.stats.pval_slope = pMin_slope;

    analysis(row).stageAvg.stageTime = stageTime;
    analysis(row).stageMax.stageTime = stageTime;
    analysis(row).stageMin.stageTime = stageTime;
    
    analysis(row).trialInfo.trialNumber = options.trialNumber{i};
    if ~isempty(options.trialTable)
        analysis(row).trialInfo.trialTable = options.trialTable(options.trialNumber{i},:);
    end
    analysis(row).options = options;
end

% Save
save(strcat(params.session.path,filesep,'analysis_',params.session.name),'analysis','-append');
disp('Finished: analysis struct created and saved');

end