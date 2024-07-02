%% Shun_analyzeSlice
% 04/11/2023
% Code partly from Bernardo

% 09/13/23
% 1. Changed file selector to uipickfiles

%% Define data path
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,colors,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;

% % Select sessions for analysis
% sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/'));
% expPath = uipickfiles('FilterSpec',parentpath,'Prompt','Select experiment folders');
% expName = erase(expPath,parentpath);

% Load excel sheet (Shun_PhysLog.xlsx)
physlog = readtable('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\wengang\Exp_withShun\Shun_PhysLog.xlsx','Sheet','Processed');

% Define result directory
resultspath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\','Select Results folder directory');

% Select experiments for analysis
parentpath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\wengang\Exp_withShun\';
expPath = uipickfiles('FilterSpec',parentpath,'Prompt','Select experiment folders');
expName = erase(expPath,parentpath);

%% Load acquisions (sweeps) for each cell

% Initalize some structures
experiments = cell(size(expPath));
epochs = cell(size(expPath));
cells = struct();
prevMouseID = ''; row = 0;

for i = 1:length(expName)
    experiments{i} = physlog(ismember(physlog.Folder,expName{i}),:);
    
    % 2. Grab corresponding recordings through epoch average files
    epochList = sortrows(struct2cell(dir(fullfile([parentpath expName{i}],['AD0_e*','p1avg.mat'])))',3);
    sweepList = cell(length(epochList),2);

    prevCellID = 0; 
    for j = 1:length(epochList)
        % Load epoch file to find individual sweep.mat
        load(fullfile([parentpath expName{i}],epochList{j,1}));
        namesplit = strsplit(epochList{j,1},{'e','p1avg'}); 
        epoch = str2double(namesplit{2});
        sweepAcq = evalin('base',['AD0_e',num2str(epoch),'p1avg.UserData.Components']);
        
        % Grab from excel sheet
        epochRow = experiments{i}(experiments{i}.Epoch == epoch, :); % Row of the corresponding epoch
        if isempty(epochRow); continue; end
        epochHeader = evalin('base',['AD0_e',num2str(epoch),'p1avg.UserData.headerString']);
        
        % Determine cellID and epoch number for that cell
        cellID = epochRow.Cell; 
        if prevCellID ~= cellID && ~strcmp(prevMouseID,epochRow.Mouse) % if we encountered a new cell
            cellEpoch = 1;
            prevCellID = cellID;
            epochs = epoch;
            v_hold = epochRow.Vhold;
            row = row + 1;
        else
            cellEpoch = cellEpoch + 1;
            epochs = [epochs, epoch];
            v_hold = [v_hold, epochRow.Vhold];
        end
        
        % Add the current cell to `cells` struct (epoch general info)
        cells(row).mouse = epochRow.Mouse{1};
        cells(row).cellID = cellID;
        cells(row).slice = epochRow.Slice;
        cells(row).v_hold = v_hold;
        cells(row).power = epochRow.Power;
        cells(row).ND_angle = epochRow.ND_angle;
        cells(row).conditions = epochRow.Conditions;
        cells(row).stim{cellEpoch} = struct();
        cells(row).sweepIdx{cellEpoch} = sweepAcq; % not working, need to find epoch number per cell
        cells(row).effects{cellEpoch} = epochRow.Effects;
        cells(row).notes{cellEpoch} = epochRow.Notes;
        cells(row).folder = epochRow.Folder;
        cells(row).epochs = epochs;
        cells(row).xscale = evalin('base',['AD0_e',num2str(epoch),'p1avg.xscale']);
        cells(row).yscale = evalin('base',['AD0_e',num2str(epoch),'p1avg.yscale']);
        cells(row).zscale = evalin('base',['AD0_e',num2str(epoch),'p1avg.zscale']);
        cells(row).outputFs = str2double(phUtil_HeaderValue(epochHeader, 'outputRate'));
        cells(row).inputFs = str2double(phUtil_HeaderValue(epochHeader, 'inputRate'));
        
        % Store Rc check params (should be the same)
        rcHeader = phUtil_HeaderValue(epochHeader, 'state.phys.internal.pulseString_RCCheck');
        cells(row).rc_check.amplitude = phUtil_parsePulsePatternString(rcHeader, 'amplitude');
        cells(row).rc_check.duration = phUtil_parsePulsePatternString(rcHeader, 'duration');
        cells(row).rc_check.offset = phUtil_parsePulsePatternString(rcHeader, 'offset');
        cells(row).rc_check.numPulses = phUtil_parsePulsePatternString(rcHeader, 'numPulses');
        cells(row).rc_check.pulseWidth = phUtil_parsePulsePatternString(rcHeader, 'pulseWidth');
        cells(row).rc_check.delay = phUtil_parsePulsePatternString(rcHeader, 'delay');
        cells(row).rc_check.isi = phUtil_parsePulsePatternString(rcHeader, 'isi');
        cells(row).rc_check.ramp = phUtil_parsePulsePatternString(rcHeader, 'ramp');
        cells(row).rc_check.patternRepeats = phUtil_parsePulsePatternString(rcHeader, 'patternRepeats');
        cells(row).rc_check.patternISI = phUtil_parsePulsePatternString(rcHeader, 'patternISI');
        
        % Initialize some temporary matrix
        sweeps = zeros(length(sweepAcq), size(evalin('base',['AD0_e',num2str(epoch),'p1avg.data']),2));
        % nx1 matrix to match with sweeps where row is each sweep
        amplitude = zeros(length(sweepAcq),1); duration = zeros(length(sweepAcq),1); 
        offset = zeros(length(sweepAcq),1); numPulses = zeros(length(sweepAcq),1); 
        pulseWidth = zeros(length(sweepAcq),1); delay = zeros(length(sweepAcq),1); 
        isi = zeros(length(sweepAcq),1); ramp = zeros(length(sweepAcq),1); 
        patternRepeats = zeros(length(sweepAcq),1); patternISI = zeros(length(sweepAcq),1); 
        vm = zeros(length(sweepAcq),1); im = zeros(length(sweepAcq),1);
        rm = zeros(length(sweepAcq),1); rs = zeros(length(sweepAcq),1);
        cm = zeros(length(sweepAcq),1);

        % Load data of individual sweeps
        for k = 1:length(sweepAcq)
            
            % Load sweep traces (.data)
            load(fullfile([parentpath expName{i}],strcat(sweepAcq{k},'.mat'))); 
            disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch)]);
            sweeps(k,:) = eval([sweepAcq{k},'.data']);
            % names{k,1} = expName{i};
            % names{k,2} = strcat(sweepAcq{k},'.mat');
            % names{k,3} = fullfile([parentpath expName{i}],strcat(sweepAcq{k},'.mat'));
            
            % Load stim patterns
            sweepHeader = eval([sweepAcq{k},'.UserData.headerString']);
            stimPattern = phUtil_HeaderValue(sweepHeader, 'state.phys.internal.pulseString_ao1');
            if length(stimPattern) < 5
                stimPattern = phUtil_HeaderValue(sweepHeader, 'state.phys.internal.pulseString_ao2');
                if length(stimPattern) < 5
                    stimPattern = phUtil_HeaderValue(sweepHeader, 'state.phys.internal.pulseString_ao0');
                end
            end
            amplitude(k) = phUtil_parsePulsePatternString(stimPattern, 'amplitude');
            duration(k) = phUtil_parsePulsePatternString(stimPattern, 'duration');
            offset(k) = phUtil_parsePulsePatternString(stimPattern, 'offset');
            numPulses(k) = phUtil_parsePulsePatternString(stimPattern, 'numPulses');
            pulseWidth(k) = phUtil_parsePulsePatternString(stimPattern, 'pulseWidth');
            delay(k) = phUtil_parsePulsePatternString(stimPattern, 'delay');
            isi(k) = phUtil_parsePulsePatternString(stimPattern, 'isi');
            ramp(k) = phUtil_parsePulsePatternString(stimPattern, 'ramp');
            patternRepeats(k) = phUtil_parsePulsePatternString(stimPattern, 'patternRepeats');
            patternISI(k) = phUtil_parsePulsePatternString(stimPattern, 'patternISI');
            
            % Load vm, im, rm, rs, cm
            vm(k) = str2double(phUtil_HeaderValue(sweepHeader, 'state.phys.cellParams.vm0'));
            im(k) = str2double(phUtil_HeaderValue(sweepHeader, 'state.phys.cellParams.im0'));
            rm(k) = str2double(phUtil_HeaderValue(sweepHeader, 'state.phys.cellParams.rm0'));
            rs(k) = str2double(phUtil_HeaderValue(sweepHeader, 'state.phys.cellParams.rs0'));
            cm(k) = str2double(phUtil_HeaderValue(sweepHeader, 'state.phys.cellParams.cm0'));    
        end
        % sweepList{epoch,1} = names;
        sweepList{epoch,2} = sweeps;
        cells(row).data{cellEpoch} = sweeps;
        
        % Store stim params for each epoch
        cells(row).stim{cellEpoch}.amplitude = amplitude;
        cells(row).stim{cellEpoch}.duration = duration;
        cells(row).stim{cellEpoch}.offset = offset;
        cells(row).stim{cellEpoch}.numPulses = numPulses;
        cells(row).stim{cellEpoch}.pulseWidth = pulseWidth;
        cells(row).stim{cellEpoch}.delay = delay;
        cells(row).stim{cellEpoch}.isi = isi;
        cells(row).stim{cellEpoch}.ramp = ramp;
        cells(row).stim{cellEpoch}.patternRepeats = patternRepeats;
        cells(row).stim{cellEpoch}.patternISI = patternISI;
        
        % Store RC check params for each epoch
        cells(row).cell_params.vm{cellEpoch} = vm;
        cells(row).cell_params.im{cellEpoch} = im;
        cells(row).cell_params.rm{cellEpoch} = rm;
        cells(row).cell_params.rs{cellEpoch} = rs;
        cells(row).cell_params.cm{cellEpoch} = cm;
        
        clearvars AD0_* sweeps names
    end
    prevMouseID = cells(end).mouse;
end
cells = sortrows_struct(cells,{'mouse','cellID'});
clearvars amplitude duration offset numPulses pulseWidth delay isi ramp patternRepeats patternISI ...
    vm im rm rs cm

%% (Not needed) Plot traces of some data

cellID = 13; cellEpoch = 2; Fs = cells(cellID).outputFs;
timeRange = 9500:12000; time = timeRange;%linspace(timeRange(1)/Fs,timeRange(end)/Fs,length(timeRange));

initializeFig(0.5,0.5);
plotSEM(time,cells(cellID).data{cellEpoch}(:,timeRange),bluePurpleRed(1,:));

%% Quality check for each cell based on epoch

cellList = 1:length(cells); % cells to analyze
peakBuffer = 2; % samples around peak loc to average
V_rev_AMPA = 0;
V_rev_GABA = -10;

for cellID = 1:length(cellList)
    disp(['Processing cell ', num2str(cellID)]);
    % Time window for baseline, stim, rc_check
    deltaT = 1/cells(cellID).xscale(2);
     
    % Sweep-specific params
    for cellEpoch = 1:length(cells(cellID).data)   
        sweepNum = size(cells(cellID).data{cellEpoch},1);
        
        % Initialize baseline properties
        I_hold_mean = zeros(sweepNum,1); I_hold_mode = zeros(sweepNum,1);
        I_hold_median = zeros(sweepNum,1); I_hold_std = zeros(sweepNum,1);
        I_hold_min = zeros(sweepNum,1); I_hold_max = zeros(sweepNum,1);
        
        % Initialize RC check properties
        rin = zeros(sweepNum,1);
        
        % Initialize PSC properties
        PSC_amplitude = zeros(sweepNum,1);
        PSC_rise_time = zeros(sweepNum,1);
        PSC_decay = zeros(sweepNum,1);
        PSC_gAMPA = zeros(sweepNum,1); PSC_gGABA = zeros(sweepNum,1);
        
        % Inialize plotting traces
        plotLength = length(950* deltaT:1200* deltaT);
        PSC_plot = zeros(sweepNum,plotLength);
        
        for k = 1:sweepNum
            disp(['Analyzing sweep ', num2str(k), ' from cellEpoch ', num2str(cellEpoch)]);
            trace = cells(cellID).data{cellEpoch}(k,:);
            
            % Define windows
            baselineWindow = (1:cells(cellID).stim{cellEpoch}.delay(k)-1) * deltaT;
            stimWindow = (cells(cellID).stim{cellEpoch}.delay(k):1200) * deltaT;
            rcWindow = cells(cellID).rc_check.delay*deltaT : length(cells(cellID).data{1});
            plotWindow = ((cells(cellID).stim{cellEpoch}.delay(k)-50)* deltaT:1200 * deltaT);
            
            % Signal processing
            % Mean subtraction
            trace_subtracted = trace - mean(trace(baselineWindow));
            Fs = 10000; % Sampling frequency  
            LP = lowpass(trace_subtracted',2000,Fs);
            % Notch filter
            d = designfilt('bandstopiir','FilterOrder',2, ...
                           'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                           'DesignMethod','butter','SampleRate',Fs);
            Notch = filtfilt(d,LP);   

            % Smooth data using sgolay filter
            yT = sgolayfilt(LP,5,27); % polynomial order of 5 and framelength of 27
            y = yT';

            % Median filter using 0.5ms window
            y = movmedian(y,6,2);

            % Subtract mean again (in SeulAh's code)
            base2 = mean(y(:,baselineWindow),2); baseM2 = repmat(base2,1,size(y,2));
            processed = y - baseM2;

            % Calculate baseline properties
            baselineTrace = processed(baselineWindow);
            I_hold_mean(k) = mean(baselineTrace);
            I_hold_mode(k) = mode(baselineTrace);
            I_hold_median(k) = median(baselineTrace);
            I_hold_std(k) = std(baselineTrace);
            I_hold_min(k) = min(baselineTrace);
            I_hold_max(k) = max(baselineTrace);
            
            % Calculate RC properties
            % rin(k) = phUtil_rcCheck_Rin();
            
            % Calculate PSC properties
            pscTrace = processed(stimWindow);
            [~,peakLoc] = max(abs(pscTrace));
            PSC_amplitude(k) = mean(pscTrace(max(peakLoc-peakBuffer,1):min(peakLoc+peakBuffer,length(pscTrace))));
            PSC_rise_time(k) = (peakLoc - stimWindow(1))/deltaT;
            %PSC_decay(k) = ;
            %ff=fit(xx(1:iPeak)', dData(1:iPeak)'-vPeak, 'exp1', 'StartPoint', [vPeak, -10] );
    %		figure; plot(ff, xx, dData-vPeak)
            %tau=-1/ff.b;
            if cells(cellID).v_hold(cellEpoch) < -50
                PSC_gAMPA(k) = (cells(cellID).v_hold(cellEpoch) - V_rev_AMPA)/-abs(PSC_amplitude(k));
            else
                PSC_gGABA(k) = (cells(cellID).v_hold(cellEpoch) - V_rev_GABA)/abs(PSC_amplitude(k));
            end
            
            % Store traces for plotting
            figTrace = processed(plotWindow);
            PSC_plot(k,:) = figTrace;
        end
        
        % Store baseline properties
        cells(cellID).cell_params.I_hold_mean{cellEpoch} = I_hold_mean;
        cells(cellID).cell_params.I_hold_mode{cellEpoch} = I_hold_mode;
        cells(cellID).cell_params.I_hold_median{cellEpoch} = I_hold_median;
        cells(cellID).cell_params.I_hold_std{cellEpoch} = I_hold_std;
        cells(cellID).cell_params.I_hold_min{cellEpoch} = I_hold_min;
        cells(cellID).cell_params.I_hold_max{cellEpoch} = I_hold_max;
        
        % Store RC properties
        
        % Store PSC properties
        cells(cellID).cell_params.PSC_amplitude{cellEpoch} = PSC_amplitude;
        cells(cellID).cell_params.PSC_rise_time{cellEpoch} = PSC_rise_time;
        cells(cellID).cell_params.PSC_decay{cellEpoch} = PSC_decay;
        cells(cellID).cell_params.PSC_gAMPA{cellEpoch} = PSC_gAMPA;
        cells(cellID).cell_params.PSC_gGABA{cellEpoch} = PSC_gGABA;
        cells(cellID).traces{cellEpoch} = PSC_plot;
    end
    % save(strcat(expPath{1},'/cells'),'cells','-v7.3');
end

save(strcat(resultspath,'/cells'),'cells','-v7.3');

%% Perform mean subtraction (baseline as 0-1000ms)

Fs = 10000;
baselineSamp = [0,1] * Fs;
baselineWindow = baselineSamp(1)+1:baselineSamp(end)-1;
processed = cell(length(sweepList),1);

for i = 1:length(sweepList)
    epochTraces = sweepList{i,2};

    if isempty(epochTraces); continue; end

    % Pre-processing
    % Calculate baseline current amplitude
    base = mean(epochTraces(:,baselineWindow),2); 
    baseM = repmat(base,1,size(epochTraces,2));
    epoch_subtracted = epochTraces - baseM;  
    
    % Low pass filter at 2kHz
    Fs = 10000; % Sampling frequency  
    LP = lowpass(epoch_subtracted',2000,Fs);
    % Notch filter
    % d = designfilt('bandstopiir','FilterOrder',2, ...
    %                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    %                'DesignMethod','butter','SampleRate',Fs);
    % Notch = filtfilt(d,LP);   
    
    % Smooth data using sgolay filter
    yT = sgolayfilt(LP,5,27); % polynomial order of 5 and framelength of 27
    y = yT';
    
    % Median filter using 0.5ms window
    y = movmedian(y,6,2);

    % figure;
    % plot(epoch_subtracted(2,:)); hold on; plot(y(2,:));

    % Subtract mean again (in SeulAh's code)
    base2 = mean(y(:,baselineWindow),2); baseM2 = repmat(base2,1,size(y,2));
    processed{i} = y - baseM2;
    
end

%% Plot PSC

animalList = unique({cells.mouse});
peakPSC = []; 
current_cutoff = 15;

for a = 1:length(animalList)
    
    cellList = find(contains({cells.mouse},animalList{a}));
    % Initialize peak arrays
    peakEPSC = nan(length(cellList),1); peakIPSC = nan(length(cellList),1);

    initializeFig(1,1); tiledlayout('flow');
    for i = 1:length(cellList)
        row = cellList(i);
        % Remove if cells only have one sweep
        if length(unique(cells(row).v_hold)) < 2; continue; end

        % Get last entry of each unique value
        [~,lastIdx] = unique(cells(row).v_hold,'last'); % sorted naturally
        firstEpoch = lastIdx(1); 
        lastEpoch = lastIdx(end);

        % Remove cells with high Rs (avg > 25)
        avgRs_exci = mean(cells(row).cell_params.rs{firstEpoch});
        avgRs_inhi = mean(cells(row).cell_params.rs{lastEpoch});
        if avgRs_exci > 50 && avgRs_inhi > 50; continue; end

        % Remove cells with <15pA response
        avgPSCAmp_exci = mean(cells(row).cell_params.PSC_amplitude{firstEpoch});
        avgPSCAmp_inhi = mean(cells(row).cell_params.PSC_amplitude{lastEpoch});
        if abs(avgPSCAmp_exci) < current_cutoff && abs(avgPSCAmp_inhi) < current_cutoff; continue; end

        % Get peak amplitude
        peakEPSC(i) = avgPSCAmp_exci; peakIPSC(i) = avgPSCAmp_inhi;

        % Plot cells
        nexttile
        plotTimeRange = 9500:12000;
        plotSEM(plotTimeRange,cells(row).traces{firstEpoch},blueWhiteRed(1,:));
        plotSEM(plotTimeRange,cells(row).traces{lastEpoch},blueWhiteRed(500,:));
        xlabel("Time (ms)"); ylabel("Current (pA)");
        title(['Cell',num2str(cells(row).cellID)]);
    end

    nexttile
    scatter(abs(peakEPSC),peakIPSC,'filled'); refline(1,0);
    xlabel("EPSC magnitude (pA)"); ylabel("IPSC magnitude (pA)");
    
    % Add to final EPSC/IPSC summary
    animalCol = ones(length(cellList),1).*a;
    peakPSC = [peakPSC; animalCol, peakEPSC, peakIPSC];
    
    saveFigures(gcf,strcat(animalList{a},"-EPSCvsIPSC"),[resultspath,'\Patch_EPSCvsIPSC']);
    close all
end

%% Plot summary plot
initializeFig(0.67,0.5); tiledlayout(1,2);
colorIdx = round(linspace(1,500,length(animalList)));
cmap = turbo(500);

nexttile;
for i = 1:length(animalList)
    animalRange = find(peakPSC(:,1) == i);
    scatter(abs(peakPSC(animalRange,2)),peakPSC(animalRange,3),'filled',...
        'MarkerFaceColor',cmap(colorIdx(i),:),'MarkerFaceAlpha',1); hold on
end
xlabel("EPSC magnitude (pA)"); ylabel("IPSC magnitude (pA)"); h1 = refline(1,0);
legend(animalList,'location','best');

nexttile;
for i = 1:length(animalList)
    animalRange = find(peakPSC(:,1) == i);
    data = rmmissing(peakPSC(animalRange,2:3));
    scatter(abs(data(:,1)),data(:,2),'filled',...
        'MarkerFaceColor',cmap(colorIdx(i),:),'MarkerFaceAlpha',1); hold on
    [p,S] = polyfit(abs(data(:,1)),data(:,2),1);
    [yfit,delta] = polyval(p,abs(data(:,1)),S);
    plot(abs(data(:,1)),yfit,'Color',[cmap(colorIdx(i),:),0.2],'LineWidth',2);
end
xlabel("EPSC magnitude (pA)"); ylabel("IPSC magnitude (pA)"); h2 = refline(1,0);
xlim([20,1400]); ylim([20,1400]); %legend(animalList); 

h1.LineWidth = 2; h2.LineWidth = 2; h1.Color = [.5,.5,.5]; h2.Color = [.5,.5,.5];
h1.LineStyle = '--'; h2.LineStyle = '--';

saveFigures(gcf,'Patch_EPSCvsIPSC_all_scatter',resultspath);
clearvars h* cmap

%% Plot distribution for all animals

