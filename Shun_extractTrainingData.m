% Shun_extractTrainingData.m

% 2025/02/07

% Extract training data to train RNN to predict photometry responses

%% Load the CSV file

clear; close all;
% sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project clamping/Recordings'));

targetFs = 200;
originalFs = 10000;
today = char(datetime('today','Format','yyyyMMdd'));

prompt = 'Enter database notes (dataset_20230326_notes.mat):';
dlgtitle = 'Save animals struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
filename = strcat('dataset_',today,'_',answer{1});

%% Loop over all select sessions and extract training data

reward_data = [];
punish_data = []; 

for s = 1:length(sessionList)

    % Load recording files
    dirsplit = split(sessionList{s},filesep);
    sessionName = dirsplit{end};
    load(strcat(sessionList{s},filesep,'data_',sessionName,'.mat'));
    load(strcat(sessionList{s},filesep,'sync_',sessionName,'.mat'));
    disp(['Extracting data from ',sessionName]);

    % Downsample photometry data
    downsampled = downsampleSignal(photometry_raw,targetFs=targetFs,originalFs=originalFs,rollingZ=false);
    photometry_ds = downsampled.dsData;
    photometry = photometry_ds; % photometry_raw

    % Find events
    lickON = find(rightLick);
    waterIdx = find(rightSolenoid);
    airpuffIdx = find(airpuff);

    waterLickIdx = nan(size(waterIdx));
    for i = 1:length(waterIdx)
        nextLick = lickON(find(lickON>=waterIdx(i),1));
        if ~isempty(nextLick); waterLickIdx(i) = nextLick; end
    end
    waterLickIdx = rmmissing(waterLickIdx);

    % Randomly selct equal amount of reward & punish events
    [rewardTraces,~] = plotTraces(waterLickIdx,[-10,10],photometry,params,...
                         signalFs=targetFs,signalSystem='NI',plot=false);
    [punishTraces,~] = plotTraces(airpuffIdx,[-10,10],photometry,params,...
                         signalFs=targetFs,signalSystem='NI',plot=false);

    % Combine multiple sessions
    reward_data = [reward_data; rewardTraces];
    punish_data = [punish_data; punishTraces];
end

% Store in .mat file
save(strcat('/Volumes/MICROSCOPE/Shun/Project clamping/Results/',filename),...
    'sessionList',"reward_data","punish_data");
disp(strcat(filename,'.mat saved'));