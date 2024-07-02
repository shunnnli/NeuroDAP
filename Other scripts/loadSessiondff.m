%% Load data_labjack.mat, data_sessionName.mat, sync_sessionName.mat, timeseries_sessionName.mat

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

%% Demodulate
signal_channel = 2;
processed = demodulateSignal(labjack.raw(signal_channel,:),...
                            targetFs=50,...
                            originalFs=labjack.samplerate,...
                            modFreq=labjack.modFreq(signal_channel));

signal_nodetrend = processed.demodData_nodetrend; % raw signal

%% Get dff traces in general

initializeFig(1,0.5);
plot(signal_nodetrend(10000:20000)/mean(signal_nodetrend));

%% Get dff from event

initializeFig(0.5,0.5);
yyaxis left
[zscore,~] = plotTraces(find(airpuff),[-0.5,3],timeSeries(2).data,params,signalFs=50,signalSystem='LJ',color=[0.124,0.434,0.73]);
yyaxis right
[dff,~] = plotTraces(find(airpuff),[-0.5,3],signal_nodetrend,params,dff=true,...
            signalFs=50,signalSystem='LJ',...
            color=[0.324,0.734,0.23],ylabel='df/f');

%% 

ds = movmean(photometry_raw,1000);
plot(ds(:)/mean(ds));
