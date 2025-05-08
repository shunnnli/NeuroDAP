% Shun_analyzeDAFFT.m

% By Shun Li and Grace Knipe, 2025

%% Load a list of sessions
clear; close all;
% addpath(genpath('/Users/graceknipe/Desktop/HMS Sabatini Lab/NeuroDAP'));
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
today = char(datetime('today','Format','yyyyMMdd'));

sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'))';

[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%% Extract FFT for all sessions

fftSummary = struct([]);

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
    [fftFreq,fftPower] = plotFFT(dLight,plot=false,print=false,...
                                 Fs=Fs);

    % (Unfinished) check whether fftFreq is the same

    fftSummary(s).animal = animalName;
    fftSummary(s).experiment = experiment;
    fftSummary(s).session = sessionName;
    fftSummary(s).power = fftPower;
    fftSummary(s).freq = fftFreq;

    disp(['Finished (',num2str(s),'/',num2str(length(sessionList)),'): ',sessionName,' analyzed']);
end

toRemove = arrayfun(@(s) all(structfun(@isempty,s)), fftSummary);
fftSummary(toRemove) = [];

%% Save fftSummary

disp('Saving fftSummary...');
rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-FFT');
resultPath = strcat(rootPath,filesep,today);
if ~isfolder(resultPath); mkdir(resultPath); end
save(fullfile(resultPath,strcat('fftSummary_',today,'.mat')),"fftSummary","-v7.3");
disp("Finished: fftSummary saved");

%% Plot fft power plot

exps = unique({fftSummary.experiment});
powers = cellfun(@(e) vertcat(fftSummary(strcmp({fftSummary.experiment},e)).power), exps, 'UniformOutput',false);
fftFreq = fftSummary(1).freq;
colorIdx = round(linspace(1,500,length(exps)));

initializeFig(0.5,1); tiledlayout(1,2);

nexttile;
for exp = 1:length(exps)
    plotSEM(fftFreq, powers{exp}, bluePurpleRed(colorIdx(exp),:),label=exps{exp});
end
xlabel('Frequency'); ylabel('Power'); legend('show');

nexttile;
for exp = 1:length(exps)
    plotSEM(fftFreq, powers{exp}, bluePurpleRed(colorIdx(exp),:),label=exps{exp});
end
xlabel('Frequency'); ylabel('Power'); legend('show');
xlim([0,3]);

%% Check whether there's any frequency differences

% number of experiments and freq‐bins
nFreq = size(powers{1},2);
colorIdx = round(linspace(1,500,length(exps)));

close all; initializeFig(0.5,1); tiledlayout(1,3);

multiCompTest = 'BH-FDR';
% multiCompTest = 'Bonferroni'; (similar results)

for exp1 = 1:length(exps)
    for exp2 = exp1+1:length(exps)
        
        expToCompare = [exp1, exp2];
        pvals = nan(1,nFreq);
        for f = 1:nFreq
            x1 = powers{exp1}(:,f);
            x2 = powers{exp2}(:,f);
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
        yyaxis right
        diff_power = mean(powers{expToCompare(1)},1) - mean(powers{expToCompare(2)},1);
        label = sprintf('%s - %s', exps{expToCompare(1)}, exps{expToCompare(2)});
        plot(fftFreq, diff_power,color=[.32 .78 .53],DisplayName=label);
        ax = gca; ax.YColor = [.32 .78 .53];
        ylabel('Power difference');

        yyaxis left
        for e = expToCompare
            % plot mean spectra for each exp
            plotSEM(fftFreq, powers{e}, bluePurpleRed(colorIdx(e),:),label=exps{e});
        end
        % mark significant freqs
        meanCells = cellfun(@(P) mean(P(:,sigIdx),1), powers, 'UniformOutput',false);
        meanMat = vertcat(meanCells{:}); % stack into an nExp×Nsig matrix
        maxVals = max(meanMat,[],1); % find the max across experiments, for each freq
        stem(fftFreq(sigIdx), maxVals*1.1, 'k','filled','DisplayName','significant');
        xlabel('Frequency (Hz)'), ylabel('Power')
        legend('Location','best')
    end
end

saveFigures(gcf,'FFT-summary',resultPath,savePNG=false,savePDF=false);

%% --- ASSUMPTION CHECK FOR T-TEST --- 
% powers{1} and powers{2} are M1×F and M2×F matrices
alpha   = 0.05;                               % significance level
nFreq   = size(powers{1},2);

% preallocate
isNormal1 = false(1,nFreq);
isNormal2 = false(1,nFreq);
equalVar  = false(1,nFreq);

for j = 1:nFreq
    x1 = powers{1}(:,j);
    x2 = powers{2}(:,j);
    
    % 1) Normality test (Lilliefors)
    isNormal1(j) = ~lillietest(x1,'Alpha',alpha);
    isNormal2(j) = ~lillietest(x2,'Alpha',alpha);
    
    % 2) Variance equality test (two‐sample F test)
    equalVar(j)  = ~vartest2(x1, x2, 'Alpha',alpha);
end

% Summarize
pctNorm1 = mean(isNormal1)*100;
pctNorm2 = mean(isNormal2)*100;
pctEqVar = mean(equalVar)*100;

fprintf('Group1 normal in %.1f%% of bins\n', pctNorm1);
fprintf('Group2 normal in %.1f%% of bins\n', pctNorm2);
fprintf('Equal variances in %.1f%% of bins\n', pctEqVar);

if pctNorm1 > 80 && pctNorm2 > 80
    disp('—> Data are approximately normal in most bins; t-test is reasonable.');
else
    disp('—> Consider a nonparametric test (e.g. ranksum or kstest2) for some bins.');
end

if pctEqVar > 80
    disp('—> Variances are equal in most bins; use ttest2(...,''Vartype'',''equal'').');
else
    disp('—> Variances differ often; use Welch''s t-test: ttest2(...,''Vartype'',''unequal'').');
end

% (Optional) Inspect a specific bin with a QQ‐plot
binToInspect = round(nFreq/2);
figure;
subplot(1,2,1), qqplot(powers{1}(:,binToInspect)), title('Group1 QQ');
subplot(1,2,2), qqplot(powers{2}(:,binToInspect)), title('Group2 QQ');

[p1,~] = lillietest(powers{1}(:,binToInspect));
[p2,~] = lillietest(powers{2}(:,binToInspect));

fprintf('Group1 normality p = %.3f\nGroup2 normality p = %.3f\n', p1, p2);


%% Test whether plotFFT is correct

% === Sample‐data FFT test script ===

% 1) Parameters
fs   = 1000;         % sampling rate (Hz)
T    = 2;            % duration (seconds)
t    = (0:1/fs:T-1/fs)';   % time vector (column)

% 2) Create a test signal:  
%    - 50 Hz sine at amplitude 1  
%    - 120 Hz sine at amplitude 0.5  
%    - added low‐level white noise
x = 1.0*sin(2*pi*50*t) + 0.5*sin(2*pi*120*t) + 0.1*randn(size(t));

% 3) Run your FFT function  
%    (replace `myFFT` with the name of your function)
%    Assume it returns two vectors: freq (Hz) and Pxx (power)
[ freq, Pxx ] = plotFFT(x, Fs=fs, xlogScale=false);

% 4) Plot the spectrum
figure;
plot(freq, Pxx, 'LineWidth',1.5)
xlim([0 200])               % zoom in on 0–200 Hz
xlabel('Frequency (Hz)')
ylabel('Power')
title('FFT of test signal')
grid on

% 5) Check that you see clear peaks at 50 Hz and 120 Hz:
%    - Peak at 50 Hz should be ~1^2/2 = 0.5 (if you’re doing a one‐sided PSD)
%    - Peak at 120 Hz should be ~(0.5)^2/2 = 0.125

% 6) (Optional) Repeat with different signals, e.g.:
%    a) Single pure tone:
%       x1 = sin(2*pi*75*t);
%    b) Two very close frequencies:
%       x2 = sin(2*pi*60*t) + 0.8*sin(2*pi*63*t);
%    c) A square wave (rich harmonics):
%       x3 = square(2*pi*30*t);

% You can wrap steps 2–4 in a loop over these variants to fully validate your FFT.
