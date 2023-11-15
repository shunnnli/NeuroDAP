function analysis = analyzeTraces(timeSeries,lick,analysisEvents,analysisLabels,analysisColors,params,options)

arguments
    timeSeries struct
    lick double
    analysisEvents cell
    analysisLabels cell
    analysisColors cell
    params

    options.timeRange double = [-15, 15]
    options.stageTime double = [-2,0;0,2]
    options.lick_binSize double = 0.1
    options.nboot double = 10000
end

disp('Ongoing: analyze traces and saved in anlaysis struct');

% Define analysis params
analysis = struct([]);
timeRange = options.timeRange;
stageTime = options.stageTime;

% Loop through all events
for i = 1:length(analysisEvents)
    %% 1. Loop through all timeSeries
    for signal = 1:size(timeSeries,2)
        disp(['     Ongoing: analyze ',analysisLabels{i},' in ',timeSeries(signal).name,' signal']);

        % Load signal of interest
        data = timeSeries(signal).data;
        finalFs = timeSeries(signal).finalFs;
        system = timeSeries(signal).system;

        % Save overall traces
        if isempty(analysisEvents{i}); continue; end
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
        stagePeak = nan(size(trace,1),size(stageTime,1));
        stagePeakFit = nan(size(stageTime,1),2);
        stagePeakFit_bs = nan(size(stageTime,1),options.nboot,2);
        pPeak_intercept = nan(size(stageTime,1),2);
        pPeak_slope = nan(size(stageTime,1),2);
        for stage = 1:size(stageTime,1)
            stageWindow = stageBin(stage,1):stageBin(stage,2);
            stagePeak(:,stage) = max(trace(:,stageWindow),[],2);

            % Fit stagePeak across session
            x = 1:size(trace,1); y = stagePeak(:,stage)';
            stagePeakFit(stage,:) = polyfit(x,y,1);

            % Test for significance (shuffled individual points, 
            % calculate slope to build distribution, and test significance)
            for sample = 1:options.nboot
                bs_data = y(randperm(length(y)));
                stagePeakFit_bs(stage,sample,:) = polyfit(x,bs_data,1);
            end
            % Calculate one-side p value
            pPeak_intercept(stage,1) = sum(stagePeakFit_bs(stage,:,end)<=stagePeakFit(stage,end))/options.nboot;
            pPeak_intercept(stage,2) = sum(stagePeakFit_bs(stage,:,end)>=stagePeakFit(stage,end))/options.nboot;
            pPeak_slope(stage,1) = sum(stagePeakFit_bs(stage,:,end-1)<=stagePeakFit(stage,end-1))/options.nboot;
            pPeak_slope(stage,2) = sum(stagePeakFit_bs(stage,:,end-1)>=stagePeakFit(stage,end-1))/options.nboot;
        end

        % Save traces and anlaysis data
        row = size(analysis,2) + 1;
        analysis(row).animal = params.session.animal;
        analysis(row).date = params.session.date;
        analysis(row).session = params.session.name;
        analysis(row).task = params.session.task;
        analysis(row).label = analysisLabels{i};
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

        analysis(row).stagePeak.data = stagePeak;
        analysis(row).stagePeak.fit = stagePeakFit;
        analysis(row).stagePeak.stats.bs = stagePeakFit_bs;
        analysis(row).stagePeak.stats.pval_intercept = pPeak_intercept;
        analysis(row).stagePeak.stats.pval_slope = pPeak_slope;

        analysis(row).stageTime = stageTime;
        analysis(row).color = analysisColors{i};
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
    stagePeak_lick = nan(size(lickRate,1),size(stageTime,1));
    stagePeakFit_lick = nan(size(stageTime,1),2);
    stagePeakFit_bs = nan(size(stageTime,1),options.nboot,2);
    pPeak_intercept = nan(size(stageTime,1),2);
    pPeak_slope = nan(size(stageTime,1),2);
    for stage = 1:size(stageTime,1)
        stageWindow = stageBin(stage,1):stageBin(stage,2);
        stagePeak_lick(:,stage) = max(lickRate(:,stageWindow),[],2);

        % Fit stageAvg across session
        x = 1:size(lickRate,1); y = stagePeak_lick(:,stage)';
        stagePeakFit_lick(stage,:) = polyfit(x,y,1);

        % Test for significance (shuffled individual points, 
        % calculate slope to build distribution, and test significance)
        for sample = 1:options.nboot
            bs_data = y(randperm(length(y)));
            stagePeakFit_bs(stage,sample,:) = polyfit(x,bs_data,1);
        end
        % Calculate one-side p value
        pPeak_intercept(stage,1) = sum(stagePeakFit_bs(stage,:,end)<=stagePeakFit(stage,end))/options.nboot;
        pPeak_intercept(stage,2) = sum(stagePeakFit_bs(stage,:,end)>=stagePeakFit(stage,end))/options.nboot;
        pPeak_slope(stage,1) = sum(stagePeakFit_bs(stage,:,end-1)<=stagePeakFit(stage,end-1))/options.nboot;
        pPeak_slope(stage,2) = sum(stagePeakFit_bs(stage,:,end-1)>=stagePeakFit(stage,end-1))/options.nboot;
    end

    % Save lick traces
    row = size(analysis,2) + 1;
    analysis(row).animal = params.session.animal;
    analysis(row).date = params.session.date;
    analysis(row).session = params.session.name;
    analysis(row).task = params.session.task;
    analysis(row).label = analysisLabels{i};
    analysis(row).name = 'Lick';
    analysis(row).system = 'Lick';
    analysis(row).data.lickRate = lickRate;
    analysis(row).data.lickEvents = lickEvents;
    analysis(row).timestamp = t_licks;
    analysis(row).timeRange = timeRange;
    analysis(row).finalFs = 1/options.lick_binSize;

    analysis(row).stageAvg.data = stageAvg_lick;
    analysis(row).stageAvg.fit = stageAvgFit_lick;
    analysis(row).stageAvg.stats.bs = stagePeakFit_bs;
    analysis(row).stageAvg.stats.pval_intercept = pAvg_intercept;
    analysis(row).stageAvg.stats.pval_slope = pAvg_slope;
    

    analysis(row).stagePeak.data = stagePeak_lick;
    analysis(row).stagePeak.fit = stagePeakFit_lick;
    analysis(row).stagePeak.stats.bs = stagePeakFit_bs;
    analysis(row).stagePeak.stats.pval_intercept = pPeak_intercept;
    analysis(row).stagePeak.stats.pval_slope = pPeak_slope;

    analysis(row).stageTime = stageTime;
    analysis(row).color = analysisColors{i};
end

% Save
save(strcat(params.session.path,filesep,'analysis_',params.session.name),'analysis','-append');
disp('Finished: analysis struct created and saved');

end