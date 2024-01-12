%%

% Shun_analyzeSessions_FFT.m

%% 
initializeFig(0.6,0.5);
nexttile; [fftFreq_NAc,P1_NAc] = plotFFT(labjack.raw(1,:)); 
title('NAc'); 
sigma_power_NAc = P1_NAc(fftFreq_NAc >= 12 & fftFreq_NAc <= 15);
peak_sigma_NAc = max(sigma_power_NAc);
avg_sigma_NAc = mean(sigma_power_NAc);

nexttile; [fftFreq_LHb,P1_LHb] = plotFFT(labjack.raw(2,:)); 
title('LHb'); 
sigma_power_LHb = P1_LHb(fftFreq_LHb >= 12 & fftFreq_LHb <= 15);
peak_sigma_LHb = max(sigma_power_LHb);
avg_sigma_LHb = mean(sigma_power_LHb);

%%

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));

sigma = struct([]);

%% Run each session
errorSessionIdx = []; errorMessage = {};

% Calculate sigma band power
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx errorMessage sigma
    
    dirsplit = strsplit(sessionList{s},filesep);
    sessionName = dirsplit{end};
    dirsplit = strsplit(sessionName,{'-','_'});
    animal = dirsplit{2};

    disp(['Calculating FFT for session ', sessionName]);

    try
        % Calculate sigma power
        load(strcat(sessionList{s},filesep,'data_labjack.mat'),'labjack');
        load(strcat(sessionList{s},filesep,'data_',sessionName,'.mat'),'photometryNI');

        sigma_power_NAc = nan; sigma_power_LHb = nan;
        peak_sigma_NAc = nan; peak_sigma_LHb = nan;
        avg_sigma_NAc = nan; avg_sigma_LHb = nan;

        for signal = 1:labjack.nSignals
            if contains(labjack.name{signal},"NAc")
                [fftFreq_NAc,P1_NAc] = plotFFT(labjack.raw(signal,:)); 
                sigma_power_NAc = P1_NAc(fftFreq_NAc >= 12 & fftFreq_NAc <= 15);
                control_low_NAc = P1_NAc(fftFreq_NAc >= 8 & fftFreq_NAc <= 11);
                control_high_NAc = P1_NAc(fftFreq_NAc >= 16 & fftFreq_NAc <= 19);
                peak_sigma_NAc = max(sigma_power_NAc);
                avg_sigma_NAc = mean(sigma_power_NAc);
                peak_control_NAc = (max(control_low_NAc)+max(control_high_NAc))/2;
                avg_control_NAc = (mean(control_low_NAc)+mean(control_high_NAc))/2;
            elseif contains(labjack.name{signal},"LHb")
                [fftFreq_LHb,P1_LHb] = plotFFT(labjack.raw(signal,:));  
                sigma_power_LHb = P1_LHb(fftFreq_LHb >= 12 & fftFreq_LHb <= 15);
                control_low_LHb = P1_LHb(fftFreq_LHb >= 8 & fftFreq_LHb <= 11);
                control_high_LHb = P1_LHb(fftFreq_LHb >= 16 & fftFreq_LHb <= 19);
                peak_sigma_LHb = max(sigma_power_LHb);
                avg_sigma_LHb = mean(sigma_power_LHb);
                peak_control_LHb = (max(control_low_LHb)+max(control_high_LHb))/2;
                avg_control_LHb = (mean(control_low_LHb)+mean(control_high_LHb))/2;
            end
        end

        % FFT on PMT signal if no LHb recorded in Labjack
        if isnan(sigma_power_LHb) && exist('photometryNI','var')
            [fftFreq_LHb,P1_LHb] = plotFFT(photometryNI);  
            sigma_power_LHb = P1_LHb(fftFreq_LHb >= 12 & fftFreq_LHb <= 15);
            control_low_LHb = P1_LHb(fftFreq_LHb >= 8 & fftFreq_LHb <= 11);
            control_high_LHb = P1_LHb(fftFreq_LHb >= 16 & fftFreq_LHb <= 19);
            peak_sigma_LHb = max(sigma_power_LHb);
            avg_sigma_LHb = mean(sigma_power_LHb);
            peak_control_LHb = (max(control_low_LHb)+max(control_high_LHb))/2;
            avg_control_LHb = (mean(control_low_LHb)+mean(control_high_LHb))/2;
        end

        % Add to corresponding row in sigma
        if isempty(sigma); animalRow = [];
        else; animalRow = find(cellfun(@(x) contains(x,animal,IgnoreCase=true), {sigma.animal})); end
        if isempty(animalRow)
            row = size(sigma,2) + 1;
            sigma(row).animal = animal;
            sigma(row).session = sessionName;
    
            sigma(row).maxNAc = peak_sigma_NAc;
            sigma(row).maxLHb = peak_sigma_LHb;
            sigma(row).avgNAc = avg_sigma_NAc;
            sigma(row).avgLHb = avg_sigma_LHb;
    
            sigma(row).maxCtrlNAc = peak_control_NAc;
            sigma(row).avgCtrlNAc = avg_control_NAc;
            sigma(row).maxCtrlLHb = peak_control_LHb;
            sigma(row).avgCtrlLHb = avg_control_LHb;
    
            sigma(row).diffMaxNAc = peak_sigma_NAc - peak_control_NAc;
            sigma(row).diffAvgNAc = avg_sigma_NAc - avg_control_NAc;
            sigma(row).diffMaxLHb = peak_sigma_LHb - peak_control_LHb;
            sigma(row).diffAvgLHb = avg_sigma_LHb - avg_control_LHb;
        else
            row = animalRow;
            sigma(row).animal = animal;
            sigma(row).session = {sigma(row).session;sessionName};
    
            sigma(row).maxNAc = [sigma(row).maxNAc;peak_sigma_NAc];
            sigma(row).maxLHb = [sigma(row).maxLHb;peak_sigma_LHb];
            sigma(row).avgNAc = [sigma(row).avgNAc;avg_sigma_NAc];
            sigma(row).avgLHb = [sigma(row).avgLHb;avg_sigma_LHb];
    
            sigma(row).maxCtrlNAc = [sigma(row).maxCtrlNAc;peak_control_NAc];
            sigma(row).avgCtrlNAc = [sigma(row).avgCtrlNAc;avg_control_NAc];
            sigma(row).maxCtrlLHb = [sigma(row).maxCtrlLHb;peak_control_LHb];
            sigma(row).avgCtrlLHb = [sigma(row).avgCtrlLHb;avg_control_LHb];
    
            sigma(row).diffMaxNAc = [sigma(row).diffMaxNAc;peak_sigma_NAc - peak_control_NAc];
            sigma(row).diffAvgNAc = [sigma(row).diffAvgNAc;avg_sigma_NAc - avg_control_NAc];
            sigma(row).diffMaxLHb = [sigma(row).diffMaxLHb;peak_sigma_LHb - peak_control_LHb];
            sigma(row).diffAvgLHb = [sigma(row).diffAvgLHb;avg_sigma_LHb - avg_control_LHb];
        end

    catch ME
        errorSessionIdx = [errorSessionIdx;s];
        msg = getReport(ME); 
        errorMessage{end+1} = msg; disp(msg);
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end

%% Save a copy of sigma
today = char(datetime('today','Format','yyyyMMdd'));
resultspath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Results');
save(strcat(resultspath,filesep,'sigma_',today),'sigma','sessionList','-v7.3');

%% Organize sigma by animal

animalList = unique({sigma.name});
sigma_byAnimal = struct([])

for a = 1:length(animalList)
    animalIdx = find(cellfun(@(x) contains(x,'SL155',IgnoreCase=true), {sigma.animal}));
    
end

%%
plot([sigma.diffMaxLHb]); hold on; plot([sigma.diffMaxNAc]);
