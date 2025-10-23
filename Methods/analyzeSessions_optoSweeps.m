function analyzeSessions_optoSweeps(sessionpath,options)

arguments
    sessionpath string
    options.outputName string %name of the output file in subfolder of the session
    
    options.laserSource string = 'clamp'
    options.durationList double = [0.1, 0.5, 1] % in sec
    options.bluePower double = [20, 30, 40]
    options.redPower double = [45, 55]

    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    
    options.plotPhotometry logical = true % Plot photometry summary plot

    options.lick_binSize double = 0.1
    options.shortTimeRange double = [-1,5]
    options.longTimeRange double = [-5,10]
end

%% Notes
% Shun_analyzeBehavior_optoPair
% Shun Li, 11/8/2022
% 02/14/2023: tidied up code, renamed to analyzeBehavior_optoPair
% 2023/07/28: packaged trial table into a function
% 2023/09/02: added camera plotting
% 2023/09/05: changed baselineIdx to selecting baseline licks 
% 2023/10/23: changed how to plot photometry signal, assume everything
% recorded in labjack

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
if ~isfield(options,'outputName'); options.outputName = dirsplit{end}; end
if ispc; projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end

% Get animal name and session date
if ~contains(options.outputName,{'-','_'})
    sessionName = dirsplit{end-1};
    dirsplit = strsplit(sessionName,{'-','_'});
else
    sessionName = options.outputName;
    dirsplit = strsplit(options.outputName,{'-','_'}); 
end
date = dirsplit{1}; animal = dirsplit{2}; sessionTask = dirsplit{3};
clear dirsplit

disp(strcat('**********',options.outputName,'**********'));
load(strcat(sessionpath,filesep,'timeseries_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'data_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'behavior_',options.outputName,'.mat'));
load(strcat(sessionpath,filesep,'sync_',options.outputName,'.mat'));

if ~isfield(params.session,'name'); params.session.name = options.outputName; end
if ~isfield(params.session,'date'); params.session.date = date; end
if ~isfield(params.session,'animal'); params.session.animal = animal; end
if ~isfield(params.session,'projectPath'); params.session.projectPath = projectPath; end

% Create analysis.mat
if ~isempty(dir(fullfile(sessionpath,"analysis_*.mat")))
    load(strcat(sessionpath,filesep,'analysis_',options.outputName,'.mat'));
else
    save(strcat(sessionpath,filesep,'analysis_',options.outputName),'sessionName','-v7.3');
    disp('Finished: analysis_.mat not found, created a new one');
end
disp(['Finished: Session ',options.outputName,' loaded']);

%% Load behaivor params

% Define baselineSystem
if ~isfield(params.session,'baselineSystem')
    params.session.baselineSystem = 'NI';
end

% Load camera signal
if params.session.withCamera
    eyeAreaIdx = find(cellfun(@(x) strcmpi(x,'eyeArea'), {timeSeries.name}));
    if ~isempty(eyeAreaIdx); eyeArea_detrend = timeSeries(eyeAreaIdx).data; end
    pupilAreaIdx = find(cellfun(@(x) strcmpi(x,'pupilArea'), {timeSeries.name}));
    if ~isempty(pupilAreaIdx); pupilArea_detrend = timeSeries(pupilAreaIdx).data; end
end

save(strcat(sessionpath,filesep,'sync_',options.outputName),'params','-append');

%% Preprocess outcome and opto data

disp('Ongoing: preprocess outcome and opto data');

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 3 10]; % in ul
% punishList = [0 0.1];
% toneList = [0 0.5 1]; % in sec
if options.round
    if ~exist('rightSolenoid_rounded','var') || options.redo
        % Round reward and tone
        rightSolenoid = rightSolenoid ./ rewardUnit;
        rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
        airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');
        
        % Save rounded data
        save(strcat(sessionpath,filesep,'timeseries_',options.outputName),'rightSolenoid_rounded','airpuff_rounded','-append');
        disp('Finished: rounding cue/outcome data');
    end
else
    rightSolenoid_rounded = rightSolenoid;
    airpuff_rounded = airpuff;
end


% Find start of opto cue
if strcmpi(options.laserSource,'clamp')
    blueLaser = blueClamp;
    redLaser = redClamp;

    % --- Extract PWM epochs (uses your single-channel extractPWMStim) ---
    [bluePatterns, ~] = extractPWMStim(blueLaser, label="blue", ...
                            duration=options.durationList, pwm=options.bluePower);
    [redPatterns,  ~] = extractPWMStim(redLaser,  label="red", ...
                            duration=options.durationList, pwm=options.redPower);

    % Save
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"bluePatterns","redPatterns","-append");
end

% Find start of lick bout
rightLickON = find(rightLick);
if ~exist('lickBout','var') || options.redo
    % Get lick bout start time (ILI < 0.5s)
    lickBout = getLickBout(rightLickON);
    save(strcat(sessionpath,filesep,'timeseries_',options.outputName),"lickBout",'-append');
end

% Combine stim&tone to form trial start
waterIdx = find(rightSolenoid_rounded);  
airpuffIdx = find(airpuff_rounded);

% Find water lick (first lick in response to water)
waterLickIdx = nan(size(waterIdx));
for i = 1:length(waterIdx)
    nextLick = rightLickON(find(rightLickON>=waterIdx(i),1));
    if ~isempty(nextLick); waterLickIdx(i) = nextLick; end
end
waterLickIdx = rmmissing(waterLickIdx);

disp('Finished: preprocess outcome and opto data');


%% Task specific params
    
% Select event idxs
toneIdx = find(leftTone);

% Select baseline idx (align to each baseline lick)
% baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{1:end-1,'BaselineLicks'}),'BaselineLicks'});
% if isempty(baselineLicks); baselineIdx = [];
% else; baselineIdx = baselineLicks(:,1); end
randomMinSample = 15*params.sync.behaviorFs;
randomMaxSample = length(params.sync.timeNI) - (15*params.sync.behaviorFs);
baselineIdx = randi([randomMinSample,randomMaxSample],100,1); % changed to rightLickON if only baseline licks


%% Create / load analysis.mat
analysisFile = fullfile(sessionpath, strcat('analysis_', options.outputName, '.mat'));

if exist(analysisFile, 'file')
    load(analysisFile); % loads sessionName and (optionally) analysis
    if ~exist('analysis','var') || options.redo
        analysis = struct('type',{},'pwm',{},'duration',{},'data',{},'timestamp',{});
    end
else
    analysis = struct('type',{},'pwm',{},'duration',{},'data',{},'timestamp',{});
    save(analysisFile,'sessionName','analysis','-v7.3');
    disp('Finished: analysis_.mat not found, created a new one');
end
disp(['Finished: Session ',options.outputName,' loaded']);

%% Plot photometry summary plots

if options.plotPhotometry && exist('timeSeries','var')
    %% Find the number of photometry channels
    photometryIdx = find(cellfun(@(x) contains(x,["NI","LJ"],"IgnoreCase",true), {timeSeries.system}));
    photometryName = cellfun(@(x) unique(x,'rows'), {timeSeries(photometryIdx).name},'UniformOutput',false);
    nSignals = length(photometryIdx);
    disp(['Finished: found ', num2str(nSignals),' photometry signals']);

    %% Plot pre-processing steps
    initializeFig(.67,.67); tiledlayout('flow');
    
    for i = 1:nSignals
        path = photometryIdx(i);
        nexttile;
        histogram(normrnd(0,1,size(timeSeries(path).data)),200); hold on
        histogram(timeSeries(path).data,200); hold on
        box off

        skew_lj = skewness(timeSeries(path).data); 
        kur_lj = kurtosis(timeSeries(path).data);
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution',timeSeries(path).name});
        
        title(timeSeries(path).name);
        subtitle(strcat("Skewness: ",num2str(skew_lj),", Kurtosis: ",num2str(kur_lj)));
    end
    
    % Save figure
    saveas(gcf,strcat(sessionpath,filesep,'Summary_photometry_distribution.png'));

    %% Loop through timeSeries struct
    for photometry = 1:nSignals

        % Load signal of interest
        path = photometryIdx(photometry);
        signal = timeSeries(path).data;
        finalFs = timeSeries(path).finalFs;
        system = timeSeries(path).system;
        name = timeSeries(path).name;
        
        %% Plot combined PSTH
        timeRange = options.shortTimeRange;
        durList = options.durationList;

        % ---------- BLUE FIGURE (columns = durations; curves = PWM freqs) ----------
        if ~isempty(bluePatterns)
            initializeFig(0.8,0.6); tiledlayout(1,3);
            for ci = 1:numel(durList)
                d = durList(ci);
                dur_mask = bluePatterns.duration_s == d;
                S = bluePatterns(dur_mask,:);
                bluePowerList = unique(S.duty);
                nBluePowers = numel(bluePowerList);

                nexttile; hold on
                alphas = linspace(0.30,1.00,nBluePowers); % low→high Hz => more solid
                for k = 1:nBluePowers
                    cur_duty = bluePowerList(k);
                    onsetTime = table2array(S(S.duty==cur_duty,'onset'));
                    c = addOpacity(bluePurpleRed(1,:), alphas(k)); % blue tint
                    [alignedMat, timestamp] = plotTraces(onsetTime, timeRange, signal, c, params, ...
                                         signalFs=finalFs, signalSystem=system);
                    analysis(end+1) = struct( ...
                        'type',    'blue', ...
                        'pwm',      cur_duty, ...
                        'duration', d, ...
                        'data',     alignedMat,...
                        'timestamp',timestamp);
                end
               
                % Black reference trace: airpuff
                [airpuffTraces,~] = plotTraces(airpuffIdx, timeRange, signal, [0 0 0], params, ...
                                   signalFs=finalFs, signalSystem=system);
                plotEvent('',d,color=bluePurpleRed(1,:));
                analysis(end+1) = struct( ...
                        'type',    'airpuff', ...
                        'pwm',      0, ...
                        'duration', 0, ...
                        'data',     airpuffTraces,...
                        'timestamp',timestamp);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(makePWMlegend(S,"Blue"));
                title(strcat('Blue (',num2str(d),' s)'));
            end
            saveas(gcf, fullfile(sessionpath, sprintf('Summary-blue-%s.png', name)));
        end

        % ---------- RED FIGURE (columns = durations; curves = PWM freqs) ----------
        if ~isempty(redPatterns)
            initializeFig(0.8,0.6); tiledlayout(1,3);
            for ci = 1:numel(durList)
                d = durList(ci);
                dur_mask = redPatterns.duration_s == d;
                S = redPatterns(dur_mask,:);
                redPowerList = unique(S.duty);
                nRedPowers = numel(redPowerList);

                nexttile; hold on
                alphas = linspace(0.30,1.00,nRedPowers); % low→high Hz => more solid
                for k = 1:nRedPowers
                    cur_duty = redPowerList(k);
                    onsetTime = table2array(S(S.duty==cur_duty,'onset'));
                    c = addOpacity(bluePurpleRed(end,:), alphas(k)); % red tint
                    [alignedMat, timestamp] = plotTraces(onsetTime, timeRange, signal, c, params, ...
                                         signalFs=finalFs, signalSystem=system);
                    analysis(end+1) = struct( ...
                        'type',    'red', ...
                        'pwm',      cur_duty, ...
                        'duration', d, ...
                        'data',     alignedMat,...
                        'timestamp',timestamp);
                end
                % Black reference trace: rewarded licks
                [waterTraces,~] = plotTraces(waterLickIdx, timeRange, signal, [0 0 0], params, ...
                                   signalFs=finalFs, signalSystem=system);
                plotEvent('',d,color=bluePurpleRed(end,:));
                analysis(end+1) = struct( ...
                        'type',    'water', ...
                        'pwm',      0, ...
                        'duration', 0, ...
                        'data',     waterTraces,...
                        'timestamp',timestamp);
                xlabel('Time (s)'); ylabel([name,' z-score']);
                legend(makePWMlegend(S,"Red"));
                title(strcat('Red (',num2str(d),' s)'));
            end
            saveas(gcf, fullfile(sessionpath, sprintf('Summary-red-%s.png', name)));
        end
    end  
end

% Append-save analysis results
analysisFile = fullfile(sessionpath, strcat('analysis_', options.outputName, '.mat'));
save(analysisFile, 'analysis', '-append');

disp('Finished: all plots and struct are plotted and saved!');
return

end