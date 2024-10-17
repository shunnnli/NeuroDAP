function plotEpochSummary(epochs,rowIdx,options)

arguments
    epochs table
    rowIdx double % epoch number to plot

    options.save logical = false
    options.saveDataPath

    options.rig string = 'Wengang'
    options.fs double = 10000

    options.dotSize double = 200
    options.passColor
    options.failColor
    options.passCtrlColor
    options.failCtrlColor
end

%% Extract epoch values

epoch = epochs{rowIdx,'Epoch'};
cell = epochs{rowIdx,'Cell'};
sweeps = epochs{rowIdx,'Raw sweeps'}{1};
protocol = epochs{rowIdx,'Protocol'}{1};
statistics = epochs{rowIdx,'Stats'}{1};
epochOptions = epochs{rowIdx,'Options'}{1};
% nSweeps = size(sweeps,1);
options.fs = epochOptions.outputFs;

% Save rig
if ~isfield(epochOptions,'rig')
    epochOptions.rig = options.rig;
else
    options.rig = epochOptions.rig;
end

% QC params
QC = epochs{rowIdx,'QC'}{1};
QCThreshold = epochOptions.QCThreshold;
% if ~isempty(fieldnames(QC)); QCIncluded = QC.included; end

% Time windows
if isfield(epochOptions,'stimDuration')
    stimDuration = epochOptions.stimDuration;
else
    stimDuration = ((protocol.numPulses(1) * protocol.isi(1))+200) * epochOptions.outputFs/1000; % 200ms recovery window after last pulse
end
if isfield(protocol,'rcCheckOnset')
    rcCheckOnset = protocol.rcCheckOnset;
    rcCheckEnd = protocol.rcCheckEnd;
else
    rcCheckOnset = 28000;
    rcCheckEnd = 30000;
end

% Get included index for scatter
included = epochs{rowIdx,'Included'}{1};
passIdx = find(included==1);
failIdx = find(included==0);

% Extract values
Rs = QC.Rs;
Rm = QC.Rm;
Cm = QC.Cm;
Ibaseline = QC.Ibaseline;
Ibaseline_std = QC.Ibaseline_std;
Verror = QC.Verror;

periMax_response = statistics.response.periMax;
periMax_control = statistics.baseline.periMax;
periMin_response = statistics.response.periMin;
periMin_control = statistics.baseline.periMin;
auc_response = statistics.response.auc;
auc_control = statistics.baseline.auc;

% Edge cases where theres no RC check
if any(isnan(Rs)); Rs = nan(length(included),1); end
if any(isnan(Rm)); Rm = nan(length(included),1); end
if any(isnan(Cm)); Cm = nan(length(included),1); end
if any(isnan(Ibaseline)); Ibaseline = nan(length(included),1); end
if any(isnan(Ibaseline_std)); Ibaseline_std = nan(length(included),1); end
if any(isnan(Verror)); Verror = nan(length(included),1); end

if any(isnan(periMax_response)); periMax_response = nan(length(included),1); end
if any(isnan(periMax_control)); periMax_control = nan(length(included),1); end
if any(isnan(periMin_response)); periMin_response = nan(length(included),1); end
if any(isnan(periMin_control)); periMin_control = nan(length(included),1); end
if any(isnan(auc_response)); auc_response = nan(length(included),1); end
if any(isnan(auc_control)); auc_control = nan(length(included),1); end

%% Initialize plotting params

dotSize = options.dotSize;

if any(included)
    if strcmp(options.rig,'Paolo'); figname = ['Cell',num2str(cell),'-Epoch',num2str(epoch),'-passed'];
    else; figname = ['Epoch',num2str(epoch),'-passed']; end
else
    if strcmp(options.rig,'Paolo'); figname = ['Cell',num2str(cell),'-Epoch',num2str(epoch),'-failed'];
    else; figname = ['Epoch',num2str(epoch),'-failed']; end
end

if ~isfield(options,'passColor'); options.passColor = [12, 173, 74]./255; end
if ~isfield(options,'failColor'); options.failColor = 1 - 0.5*(1-options.passColor); end
if ~isfield(options,'passCtrlColor'); options.passCtrlColor = [.7 .7 .7]; end
if ~isfield(options,'failCtrlColor'); options.failCtrlColor = [.9 .9 .9]; end

%% Plot figure

initializeFig(1,1); tiledlayout(3,9); 

% Plot good acquisitions (whole trace)
nexttile(1,[1,3]);
plotWindow = 1:size(sweeps,2);
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==1,plotWindow),options.passColor,plotIndividual=true,IndividualColor=options.passCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': whole trace (included sweeps)']);

% Plot good acquisitions (stim period)
nexttile(4,[1,3]);
timeRangeStartSample = protocol.stimOnset(1) - 200;
timeRangeEndSample = timeRangeStartSample + stimDuration;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==1,plotWindow),options.passColor,plotIndividual=true,IndividualColor=options.passCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': stim trace (included sweeps)']);

% Plot good acquisitions (RC period)
nexttile(7,[1,3]);
timeRangeStartSample = rcCheckOnset - 500;
timeRangeEndSample = rcCheckEnd;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==1,plotWindow),options.passColor,plotIndividual=true,IndividualColor=options.passCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': RC trace (included sweeps)']);

% Plot bad acquisitions (whole trace)
nexttile(10,[1,3]);
plotWindow = 1:size(sweeps,2);
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==0,plotWindow),options.failColor,plotIndividual=true,IndividualColor=options.failCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': whole trace (removed sweeps)']);

% Plot bad acquisitions (stim period)
nexttile(13,[1,3]);
timeRangeStartSample = protocol.stimOnset(1) - 200;
timeRangeEndSample = timeRangeStartSample + stimDuration;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==0,plotWindow),options.failColor,plotIndividual=true,IndividualColor=options.failCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': stim trace (removed sweeps)']);

% Plot bad acquisitions (RC period)
nexttile(16,[1,3]);
timeRangeStartSample = rcCheckOnset - 500;
timeRangeEndSample = rcCheckEnd;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = plotWindow ./ (options.fs/1000);
plotSEM(timeRangeInms,sweeps(included==0,plotWindow),options.failColor,plotIndividual=true,IndividualColor=options.failCtrlColor);
xlabel('Time (ms)'); xlim([timeRangeInms(1),timeRangeInms(end)]);
ylabel('I (pA)');
title(['Epoch ',num2str(epoch),': RC trace (removed sweeps)']);

% Plot Rs vs acq
nexttile(19);
yline(QCThreshold.Rs,'--','Threshold'); hold on;
scatter(failIdx,Rs(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,Rs(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('Rs (M\Omega)');
title('Rs');

% Plot Rm vs acq
nexttile(20);
scatter(failIdx,Rm(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,Rm(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('Rm (M\Omega)');
title('Rm');

% Plot Cm vs acq
nexttile(21);
scatter(failIdx,Cm(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,Cm(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('Cm (uF)');
title('Cm');

% Plot Ibaseline vs acq
nexttile(22);
yline(QCThreshold.Ibaseline,'--','Threshold'); hold on;
scatter(failIdx,Ibaseline(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,Ibaseline(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('I (pA)');
title('Baseline current');

% Plot Ibaseline_std vs acq
nexttile(23);
yline(QCThreshold.Ibaseline_std,'--','Threshold'); hold on;
scatter(failIdx,Ibaseline_std(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,Ibaseline_std(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('Standard deviation');
title('SD of baseline current');

% Plot Verror vs acq
nexttile(24);
yline(QCThreshold.Verror,'--','Threshold'); hold on;
scatter(failIdx,abs(Verror(failIdx)),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,abs(Verror(passIdx)),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('|Voltage error| (mV)');
title('Voltage error');

% Plot max vs acq (stim & baseline)
nexttile(25);
scatter(failIdx,periMax_control(failIdx),dotSize,options.failCtrlColor,'filled'); hold on;
scatter(passIdx,periMax_control(passIdx),dotSize,options.passCtrlColor,'filled'); hold on;
scatter(failIdx,periMax_response(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,periMax_response(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('I_{max}  (pA)');
title('Max response');

% Plot min vs acq (stim & baseline)
nexttile(26);
scatter(failIdx,periMin_control(failIdx),dotSize,options.failCtrlColor,'filled'); hold on;
scatter(passIdx,periMin_control(passIdx),dotSize,options.passCtrlColor,'filled'); hold on;
scatter(failIdx,periMin_response(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,periMin_response(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('I_{min}  (pA)');
title('Min response');

% Plot charge vs acq (stim & baseline)
nexttile(27);
scatter(failIdx,auc_control(failIdx),dotSize,options.failCtrlColor,'filled'); hold on;
scatter(passIdx,auc_control(passIdx),dotSize,options.passCtrlColor,'filled'); hold on;
scatter(failIdx,auc_response(failIdx),dotSize,options.failColor,'filled'); hold on;
scatter(passIdx,auc_response(passIdx),dotSize,options.passColor,'filled'); hold on;
xlabel('Acquisition');
ylabel('Charge (pC)');
title('Charge');

if options.save
    % Find the newest results folder
    if ~isfield(options,'saveDataPath')
        sessionPath = epochs{1,'Session'};
        resultsFolders = sortrows(struct2cell(dir(fullfile(sessionPath,"Epochs-*")))',3);
        resultFolder = resultsFolders{end,1};
        options.saveDataPath = fullfile(sessionPath,resultFolder);
    end

    % Save figures
    saveFigures(gcf,figname,options.saveDataPath);
end

end