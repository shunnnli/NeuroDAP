% run_salt
% MGC 4/12/2022

paths = struct;
paths.npx_data = 'E:\npx_data\with_cam\';
paths.save = 'I:\My Drive\UchidaLab\DA_independence\salt_output\';

opt = struct;
opt.session = 'MC25_20211021.mat';

% trial events
opt.las_start = 2; % seconds after trial start
opt.iti_start = 3; % seconds after trial start
opt.iti_win = [3 5]; % seconds within iti to use for SALT baseline


% for defining good cells
opt.pres_ratio = 0.9;
opt.isi_viol = 0.5;
opt.amp_cutoff = 0.1;

% SALT params:
opt.dt = 0.001;
opt.wn = 0.02;

opt.make_psth = false;

%% load data
load(fullfile(paths.npx_data,opt.session));

%% get good cells

pass_qc = sp.cgs'==2 & ...
    sp.metrics.presence_ratio > opt.pres_ratio & ...
    sp.metrics.isi_viol < opt.isi_viol & ...
    sp.metrics.amplitude_cutoff < opt.amp_cutoff;

cellid = sp.cids(pass_qc);

%% get trigger times

trig_las = 1000 * (SessionData.TrialStartTimestamp + opt.las_start);
trig_iti = 1000 * (SessionData.TrialStartTimestamp + opt.iti_start);
trial_type = reorder_trial_types(SessionData.TrialTypes,exp_params);
nTrialTypes = numel(unique(trial_type));

nPulse = SessionData.TrialSettings(1).NumLaserPulse;
PulseDur = SessionData.TrialSettings(1).LaserPulseDuration;
PulseFreq = SessionData.TrialSettings(1).LaserPulseFrequency;

PulseDur_ms = 1000 * PulseDur;
TotalDur_ms = 1000 * (1/PulseFreq);

assert(PulseDur==opt.wn); % Assumes test window for SALT is the duration of the laser pulse


%% iterate over cells
salt_p = nan(numel(cellid),nPulse,nTrialTypes);

parfor cIdx = 1:numel(cellid)
    fprintf('Cell %d/%d: %d\n',cIdx,numel(cellid),cellid(cIdx));
    
    spiket_ms = 1000 * sp.st(sp.clu==cellid(cIdx));

    
    for ttIdx = 1:nTrialTypes

        % get discretized spike rasters
        pre_ms_test = 10;
        [~,psth_baseline] = plot_timecourse('timestamp',spiket_ms,trig_iti(trial_type==ttIdx)',...
            1000*opt.iti_win(1),1000*opt.iti_win(2),[],'win_len',1,'resample_bin',1,'plot_type','none');
        [~,psth_test] = plot_timecourse('timestamp',spiket_ms,trig_las(trial_type==ttIdx)',...
            -pre_ms_test,nPulse*TotalDur_ms,[],'win_len',1,'resample_bin',1,'plot_type','none');

        spt_baseline = psth_baseline.rate_rsp/1000;
        spt_test = psth_test.rate_rsp/1000;

        spt_baseline = spt_baseline(:,1:(1000*diff(opt.iti_win)));

        % run salt for each laser pulse
        for pIdx = 1:nPulse
            winIdx = pre_ms_test + (1:PulseDur_ms) + TotalDur_ms*(pIdx-1);
            salt_p(cIdx,pIdx,ttIdx) = salt(spt_baseline,spt_test(:,winIdx),opt.dt,opt.wn);
        end   
    end

    if opt.make_psth
        figure;
        plot_timecourse('timestamp',spiket_ms,trig_las',-TotalDur_ms,nPulse*TotalDur_ms,trial_type','win_len',1,'resample_bin',1)
%         title(sprintf('c%d: p1=%0.3f, p2=%0.3f, p3=%0.3f, p4=%0.3f', ...
%             cellid(cIdx), salt_p(cIdx,1,1), salt_p(cIdx,1,2), salt_p(cIdx,1,3), salt_p(cIdx,1,4)));
    end

end


%% save results
save(fullfile(paths.save,opt.session),'salt_p','cellid','opt');