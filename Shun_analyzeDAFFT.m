% Shun_analyzeDAFFT.m

% By Shun Li and Grace Knipe, 2025

%% Load a list of sessions
clear; close all;
% addpath(genpath('/Users/graceknipe/Desktop/HMS Sabatini Lab/NeuroDAP'));
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'))';

[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%% Extract FFT for all sessions

fftSummary = struct([]);
fftFreq = linspace(0,25,2500);

for s = 1:length(sessionList)
    sessionpath = sessionList{s};  
    
    % Extract session name from path
    dirsplit = strsplit(sessionpath, filesep); 
    sessionName = dirsplit{end};

    % Extract animal name
    animalName = strsplit(sessionName,'-'); 
    animalName = animalName{end-1};

    % Extract experiment
    experiment = dirsplit{end-1};

    % Load timesereis.mat file
    load(fullfile(sessionpath, ['timeseries_', sessionName, '.mat']));
   
    % extract dLight signal
    row = find(strcmp({timeSeries.name},'dLight'));
    if isempty(row); continue; end
    dLight = timeSeries(row).data;
    Fs = timeSeries(row).finalFs;

    % Calculate FFT on the whole DA recording
    [fftFreq_raw,fftPower_raw] = plotFFT(dLight,Fs=Fs,timeToEstimateCarrier=100,...
                                         plot=false,print=false);
    fftPower = interp1(fftFreq_raw,fftPower_raw,fftFreq,'pchip',NaN);

    fftSummary(s).animal = animalName;
    fftSummary(s).experiment = experiment;
    fftSummary(s).session = sessionName;
    fftSummary(s).power = fftPower;
    fftSummary(s).freq = fftFreq;

    % Autoregression model
    data = dLight(2:end)';           % ensure column
    T    = numel(data);
    orders = 1:2;             % model orders to compare
    
    for p = orders
        [a,e]   = aryule(data, p);        % a = [1, -φ₁, …, -φₚ], e = noise var
        phi{p}  = -a(2:end);              % extract φ₁…φₚ
        ll      = -T/2*(log(2*pi*e) + 1); % Gaussian log-likelihood
        AIC(p)  = -2*ll + 2*p;
        BIC(p)  = -2*ll + p*log(T);
    end

    ar1.phi = phi{1}; ar1.aic = AIC(1); ar1.bic = BIC(1);
    ar2.phi = phi{2}; ar2.aic = AIC(2); ar2.bic = BIC(2);

    fftSummary(s).AR1 = ar1;
    fftSummary(s).AR2 = ar2;

    disp(['Finished (',num2str(s),'/',num2str(length(sessionList)),'): ',sessionName,' analyzed']);
end

toRemove = arrayfun(@(s) all(structfun(@isempty,s)), fftSummary);
fftSummary(toRemove) = [];

%% Save fftSummary

disp('Saving fftSummary...');
rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-FFT');
today = char(datetime('today','Format','yyyyMMdd'));
resultPath = strcat(rootPath,filesep,today);
if ~isfolder(resultPath); mkdir(resultPath); end
save(fullfile(resultPath,strcat('fftSummary_',today,'.mat')),"fftSummary","-v7.3");
disp("Finished: fftSummary saved");

%% Check whether there's any frequency differences

% Extract data
exps = unique({fftSummary.experiment});
fftFreq = linspace(0,25,2500);
powers = cellfun(@(e) vertcat(fftSummary(strcmp({fftSummary.experiment},e)).power), exps, 'UniformOutput',false);
colorIdx = flip(round(linspace(1,500,length(exps))));

% Plot settings
logscale = true;
xlimit = [5,25];
smoothing = false;
normalize = true;
multiCompTest = 'BH-FDR';
% multiCompTest = 'Bonferroni'; (similar results)

% Process powers
if normalize
    powers_norm = cell(size(powers));
    for e = 1:numel(powers)
        P = powers{e};                 % M×F
        normFactor = sum(P,2);         % M×1 vector of totals
        powers_norm{e} = P ./ normFactor;   % broadcasts so each row sums to 1
    end
else
    powers_norm = powers;
end
if smoothing
    windowSize = 1;
    powers_smoothed = cellfun(@(P) movmean(P, windowSize, 2), ...
                          powers_norm, 'UniformOutput', false);
else
    powers_smoothed = powers_norm;
end


% Plot
nFreq = size(powers_smoothed{1},2);
close all; initializeFig(0.5,1); tiledlayout(1,3);

for exp1 = 1:length(exps)
    for exp2 = exp1+1:length(exps)
        
        expToCompare = [exp1, exp2];
        pvals = nan(1,nFreq);
        for f = 1:nFreq
            x1 = powers_smoothed{exp1}(:,f);
            x2 = powers_smoothed{exp2}(:,f);
            [~,pvals(f)] = kstest2(x1,x2);
        end
        
        if contains(multiCompTest,'FDR')
            % Benjamini–Hochberg FDR (no toolbox)
            alpha = 0.05;
            m = numel(pvals);
            [p_sorted, sortIdx] = sort(pvals);
            BH_thresh = (1:m)/m * alpha;
            rejections = p_sorted <= BH_thresh;
            max_i = find(rejections,1,'last');
            
            sigIdx = false(1,m);
            if ~isempty(max_i)
                sigIdx(sortIdx(1:max_i)) = true;
            end
            sigFreqs = fftFreq(sigIdx);
        elseif contains(multiCompTest,'Bon')
            % assume pvals is 1×F vector of raw p‐values, freq is 1×F frequency vector
            alpha    = 0.05;         % your family‐wise α
            m        = numel(pvals); % number of tests (frequency bins)
            
            % --- Bonferroni correction ---
            bonThresh = alpha/m;     % each test must be below α/m
            sigIdx    = pvals < bonThresh;  % logical mask of significant bins
            sigFreqs  = fftFreq(sigIdx);       % those frequencies
        end
        
        nexttile;
        % yyaxis right
        % diff_power = mean(powers_smoothed{expToCompare(1)},1) - mean(powers_smoothed{expToCompare(2)},1);
        % label = sprintf('%s - %s', exps{expToCompare(1)}, exps{expToCompare(2)});
        % plot(fftFreq, diff_power,color=[.32 .78 .53],DisplayName=label); hold on
        % % scatter(fftFreq(sigIdx), diff_power(sigIdx)*1.1, 50, 'k','filled',...
        % %         MarkerFaceAlpha=0.5, DisplayName='significant');
        % ax = gca; ax.YColor = [.32 .78 .53]; 
        % ylabel('Power difference');
        % yyaxis left

        for e = expToCompare
            % plot mean spectra for each exp
            plotSEM(fftFreq, powers_smoothed{e}, bluePurpleRed(colorIdx(e),:),label=exps{e});
        end
        % mark significant freqs
        meanCells = cellfun(@(P) mean(P(:,sigIdx),1), powers_smoothed, 'UniformOutput',false);
        meanMat = vertcat(meanCells{:}); % stack into an nExp×Nsig matrix
        maxVals = max(meanMat,[],1); % find the max across experiments, for each freq
        scatter(fftFreq(sigIdx), maxVals*1.2, 50, 'k','filled',...
                MarkerFaceAlpha=0.5, DisplayName='significant');
        xlabel('Frequency (Hz)'); ylabel('Power');
        ax = gca; ax.YColor = [0 0 0];
        legend('Location','best')
        
        xlim(xlimit);
        if logscale; set(gca, 'YScale', 'log'); end
    end
end

saveFigures(gcf,'FFT-summary',resultPath,savePNG=false,savePDF=false);

%% Autoregression model: compare second-lag coeff

exps = unique({fftSummary.experiment});
ar_results = struct([]);

for i = 1:numel(exps)
    mask              = strcmp({fftSummary.experiment}, exps{i});
    S                 = fftSummary(mask);
    ar1               = [S.AR1];        % 1×n struct array
    ar2               = [S.AR2];        % 1×n struct array

    ar_results(i).experiment = exps{i};
    ar_results(i).phi1  = cell2mat({ar1.phi})';         % [φ1 φ1 …] (1×n)
    ar_results(i).aic1  = cell2mat({ar1.aic})';
    ar_results(i).bic1  = cell2mat({ar1.bic})';

    % each ar2.phi is a 1×2 vector, so we reshape into n×2
    phi2_mat  = reshape([ar2.phi], numel(ar2(1).phi), [])';
    ar_results(i).phi2 = phi2_mat;         % n×2

    ar_results(i).aic2 = cell2mat({ar2.aic})';
    ar_results(i).bic2 = cell2mat({ar2.bic})';
end

initializeFig(0.5,1); tiledlayout(1,3);

% Compare phi2
nexttile;
taCasp3 = ar_results(1).phi2(:,2);
paAIP2 = ar_results(2).phi2(:,2);
iGluSnFR = ar_results(3).phi2(:,2);
plotScatterBar(1,taCasp3,color=bluePurpleRed(colorIdx(1),:));
plotScatterBar(2,iGluSnFR,color=bluePurpleRed(colorIdx(3),:));
plotScatterBar(3,paAIP2,color=bluePurpleRed(colorIdx(2),:));
plotStats(taCasp3,iGluSnFR,[1 2],testType='kstest');
plotStats(iGluSnFR,paAIP2,[2 3],testType='kstest');
plotStats(taCasp3,paAIP2,[1 3],testType='kstest');
xticks([1 2 3]); xticklabels({'taCasp3','iGluSnFR','paAIP2'});
ylabel('Second lag coefficient');

% Compare aic2
nexttile;
taCasp3 = ar_results(1).aic2 - ar_results(1).aic1;
paAIP2 = ar_results(2).aic2 - ar_results(2).aic1;
iGluSnFR = ar_results(3).aic2 - ar_results(3).aic1;
plotScatterBar(1,taCasp3,color=bluePurpleRed(colorIdx(1),:));
plotScatterBar(2,iGluSnFR,color=bluePurpleRed(colorIdx(3),:));
plotScatterBar(3,paAIP2,color=bluePurpleRed(colorIdx(2),:));
plotStats(taCasp3,iGluSnFR,[1 2],testType='kstest');
plotStats(iGluSnFR,paAIP2,[2 3],testType='kstest');
plotStats(taCasp3,paAIP2,[1 3],testType='kstest');
xticks([1 2 3]); xticklabels({'taCasp3','iGluSnFR','paAIP2'});
ylabel('AIC2 - AIC1');

% Compare BIC
nexttile;
taCasp3 = ar_results(1).bic2 - ar_results(1).bic1;
paAIP2 = ar_results(2).bic2 - ar_results(2).bic1;
iGluSnFR = ar_results(3).bic2 - ar_results(3).bic1;
plotScatterBar(1,taCasp3,color=bluePurpleRed(colorIdx(1),:));
plotScatterBar(2,iGluSnFR,color=bluePurpleRed(colorIdx(3),:));
plotScatterBar(3,paAIP2,color=bluePurpleRed(colorIdx(2),:));
plotStats(taCasp3,iGluSnFR,[1 2],testType='kstest');
plotStats(iGluSnFR,paAIP2,[2 3],testType='kstest');
plotStats(taCasp3,paAIP2,[1 3],testType='kstest');
xticks([1 2 3]); xticklabels({'taCasp3','iGluSnFR','paAIP2'});
ylabel('BIC2 - BIC1');

%% 

TD = [353.12 193.33 363.6 329.55 347.55 227.64 409.2 281.17 602.36 338.43];
PD = [229.27 219.9 298.58 377.30 355.2 206.07 321.79 437.2 474.97 354.12];

plotScatterBar(1,TD);
plotScatterBar(2,PD);
plotStats(TD,PD,[1 2]);


%% Load event align traces

% keep = ismember({animals.name}, {'iGluSnFR','dLight'});
% animalsFFT = animals(keep);
% combined_FFT = [combined_FFT, animalsFFT];

% [combined_FFT.experiment] = deal([]);
% [combined_FFT(1:215).experiment] = deal('iGluSnFR');
% [combined_FFT(216:338).experiment] = deal('taCasp3');
% [combined_FFT(339:end).experiment] = deal('paAIP2');
% 
% animals_iGluSnFR = combined_FFT(1:215);
% animals_taCasp3 = combined_FFT(216:338);
% animals_paAIP2 = combined_FFT(339:end);

% Save animals struct
% prompt = 'Enter database notes (animals_20230326_notes.mat):';
% dlgtitle = 'Save animals struct'; fieldsize = [1 45]; definput = {''};
% answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
% today = char(datetime('today','Format','yyyyMMdd'));
% filename = strcat('animals_fft_',today,'_',answer{1});
% resultspath = strcat('/Volumes/MICROSCOPE/Shun/Project valence/Results/Summary-FFT');
% 
% % Save animals.mat
% if ~isempty(answer)
%     disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
%     save(strcat(resultspath,filesep,filename),'animals_iGluSnFR','animals_taCasp3','animals_paAIP2','-v7.3');
%     disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
% end

% Load files
filepath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-FFT'))';
load(filepath{1});

% Set color
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
iGluSnFRColor = bluePurpleRed(1,:);
taCasp3Color = bluePurpleRed(250,:);
paAIP2Color = bluePurpleRed(500,:);

%% Compare event align traces across experiments

% Water
initializeFig(0.3,0.5); 
event = 'reward';
task = 'All';
signal = 'dLight';
totalTrialRange = [1,10];

close all; initializeFig(0.5,0.5); tiledlayout('flow');
combinedTraces_iGluSnFR = combineTraces(animals_iGluSnFR,timeRange=[-0.5,3],...
                            eventRange=event,taskRange=task,signalRange=signal,...
                            totalTrialRange=totalTrialRange,combineStats=false);
if contains(signal,'dLight')
    combinedTraces_taCasp3 = combineTraces(animals_taCasp3,timeRange=[-0.5,3],...
                                eventRange=event,taskRange=task,signalRange=signal,...
                                totalTrialRange=totalTrialRange,combineStats=false);
    combinedTraces_paAIP2 = combineTraces(animals_paAIP2,timeRange=[-0.5,3],...
                                eventRange=event,taskRange=task,signalRange=signal,...
                                totalTrialRange=totalTrialRange,combineStats=false);
end

timestamp = combinedTraces_iGluSnFR.timestamp;
data_iGluSnFR = combinedTraces_iGluSnFR.data{1};
if contains(signal,'dLight')
    data_taCasp3 = combinedTraces_taCasp3.data{1};
    data_paAIP2 = combinedTraces_paAIP2.data{1};
end

nexttile;
plotTraces(data_iGluSnFR,timestamp,color=iGluSnFRColor);
if contains(signal,'dLight')
    plotTraces(data_taCasp3,timestamp,color=taCasp3Color);
    plotTraces(data_paAIP2,timestamp,color=paAIP2Color);
end
plotEvent('Water',0,color=bluePurpleRed(1,:))
xlabel('Time (s)'); ylabel('z-score');

% FFT on event aligned trace
nexttile;
[fftFreq,fft_iGluSnFR] = plotFFT(data_iGluSnFR,color=iGluSnFRColor,Fs=50,print=false);
if contains(signal,'dLight')
    [~,fft_taCasp3] = plotFFT(data_taCasp3,color=taCasp3Color,Fs=50,print=false);
    [~,fft_paAIP2] = plotFFT(data_paAIP2,color=paAIP2Color,Fs=50,print=false);
end
xlabel('Frequency (Hz)'); ylabel('Power'); xlim([5,25]);

%% ACF analysis (sensor too slow to see anything)

exps    = {'iGluSnFR','taCasp3','paAIP2'};
allData  = {data_iGluSnFR, data_taCasp3, data_paAIP2};
colors = {iGluSnFRColor; taCasp3Color; paAIP2Color};
maxLag   = 200;   % how many lags to show (in samples)
Fs = 50;

% Preallocate
nCond  = numel(exps);
allACF = cell(nCond,1);
lags   = (0:maxLag)/Fs;

% Compute trial‐wise ACFs and then average
for c = 1:nCond
    X      = allData{c};           % nTrials × T
    [nTrials,T] = size(X);
    acfs   = zeros(nTrials, maxLag+1);
    
    for i = 1:nTrials
        % full autocorr from -maxLag..+maxLag, normalized
        r = xcorr(X(i,:), maxLag, 'coeff');
        % keep only non-negative lags (center is at maxLag+1)
        acfs(i,:) = r(maxLag+1 : end);
    end
    
    allACF{c} = acfs;
end

% --- Plotting ---
initializeFig(.5,.5);
for c = 1:nCond
    plotSEM(lags,allACF{c},colors{c});
end
xlabel('Lag (s)');
ylabel('Autocorrelation');
legend(exps,'Location','northeast');


%% Cross correlation (not useful)

close all; initializeFig(0.7,0.5); tiledlayout('flow');
LineWidth = 4; 
% rng(0, 'twister');

nexttile;
c = crossCorr2D(data_iGluSnFR,data_taCasp3,average=true);
imagesc(timestamp, timestamp, c); colorbar;
xline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
yline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
plot(timestamp, timestamp,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
xlabel('taCasp3: time (s)');
ylabel('iGluSnFR: time (s)');
title('iGluSnFR vs taCasp3');

nexttile;
c = crossCorr2D(data_iGluSnFR,data_paAIP2,average=true);
imagesc(timestamp, timestamp, c); colorbar;
xline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
yline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
plot(timestamp, timestamp,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
xlabel('paAIP2: time (s)');
ylabel('iGluSnFR: time (s)');
title('iGluSnFR vs paAIP2');

nexttile;
c = crossCorr2D(data_taCasp3,data_paAIP2,average=true);
imagesc(timestamp, timestamp, c); colorbar;
xline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
yline(0,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
plot(timestamp, timestamp,LineWidth=LineWidth,Color='w',LineStyle='--'); hold on;
xlabel('paAIP2: time (s)');
ylabel('taCasp3: time (s)');
title('taCasp3 vs paAIP2');
