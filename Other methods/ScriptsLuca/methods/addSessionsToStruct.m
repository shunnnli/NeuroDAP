function [targetStruct] = addSessionsToStruct(addlist, targetStruct, task, options)

arguments
    addlist cell % sessions to add
    targetStruct struct
    task string % baseline, baseline->reward, reward->punish, punish->reward
    options.newStruct logical = false
    options.maxMissesAllowed = 3;
    options.performanceEvalLength = 10;
end

% Examine how many sessions are already in the targetStruct
if options.newStruct
    structlen = 0;
else
    structlen = length(targetStruct);
end

% Loading all data into cell
if strcmp(task,'baseline')
    for i = 1:length(addlist)
        % Get session idx
        sessionID = structlen + i;

        % Get session information
        dirsplit = strsplit(addlist{i},filesep); 
        sessionName = dirsplit{end};
        namesplit = strsplit(sessionName,'-');
        sessionAnimal = namesplit{2};

        % Store session information in struct
        targetStruct(sessionID).mouse = sessionAnimal;
        targetStruct(sessionID).session = sessionName;
        targetStruct(sessionID).type = task;

        % Load recording data
        % Some session might not have LJ or NI photometry, initialize related
        % variable to 0 to avoid error
        rollingGreenLP = 0; timePhotometry = 0; photometryNI = 0;
        trials = 0; blueLaser = 0; airpuff_rounded = 0; firstPulse = 0;
        load(strcat(addlist{i},'\','sync_',sessionName,'.mat'),...
            'trials','params','rightLick','blueLaser',...
            'timeNI','timePhotometry','rollingGreenLP','photometryNI',...
            'rightSolenoid','leftTone','airpuff_rounded','airpuff',...
            'firstPulse','redLaser');
        if airpuff_rounded == 0; airpuff_rounded = airpuff; end
        if firstPulse == 0; firstPulse = find(redLaser); end

        % For old sessions
        if ~isfield(params.sync,'timeNI'); params.sync.timeNI = timeNI; end
        if ~isfield(params.sync,'timePhotometry'); params.sync.timePhotometry = timePhotometry; end
        if ~isfield(params.sync,'behaviorFs'); params.sync.behaviorFs = 10000; end

        % Store session recordings in struct
        targetStruct(sessionID).params = params;
        targetStruct(sessionID).trials = trials;
        targetStruct(sessionID).licks = find(rightLick==1);
        targetStruct(sessionID).water = find(rightSolenoid);
        targetStruct(sessionID).airpuff = find(airpuff_rounded);
        targetStruct(sessionID).tone = find(leftTone);
        targetStruct(sessionID).stim = firstPulse;
        targetStruct(sessionID).randomShutter = round(rand([500,1])*length(blueLaser));
        targetStruct(sessionID).photometryLJ = rollingGreenLP;
        targetStruct(sessionID).photometryNI = photometryNI;

        % Initialize analysis params (window to analyze)
        % Behavior analysis window
        targetStruct(sessionID).analysis.behavior.water = 1:length(targetStruct(sessionID).water);
        targetStruct(sessionID).analysis.behavior.airpuff = 1:length(targetStruct(sessionID).airpuff);
        targetStruct(sessionID).analysis.behavior.stim = 1:length(targetStruct(sessionID).stim);
        targetStruct(sessionID).analysis.behavior.tone = 1:length(targetStruct(sessionID).tone);
        % DA analysis window
        targetStruct(sessionID).analysis.DA.water = 1:length(targetStruct(sessionID).water);
        targetStruct(sessionID).analysis.DA.airpuff = 1:length(targetStruct(sessionID).airpuff);
        targetStruct(sessionID).analysis.DA.stim = 1:length(targetStruct(sessionID).stim);
        targetStruct(sessionID).analysis.DA.tone = 1:length(targetStruct(sessionID).tone);
        % LHb analysis window
        targetStruct(sessionID).analysis.LHb.water = 1:length(targetStruct(sessionID).water);
        targetStruct(sessionID).analysis.LHb.airpuff = 1:length(targetStruct(sessionID).airpuff);
        targetStruct(sessionID).analysis.LHb.stim = 1:length(targetStruct(sessionID).stim);
        targetStruct(sessionID).analysis.LHb.tone = 1:length(targetStruct(sessionID).tone);

        clear dirsplit namesplit params rollingGreenLP photometryNI ...
            trials rightLick airpuff_rounded timeNI timePhotometry leftTone ...
            rightSolenoid firstPulse airpuff redLaser blueLaser
        disp(['Session ',sessionName,' loaded']);
    end
    
else
    for i = 1:length(addlist)
        % Get session idx
        sessionID = structlen + i;
        
        % Get session information
        dirsplit = strsplit(addlist{i},filesep); 
        sessionName = dirsplit{end};
        namesplit = strsplit(sessionName,'-');
        sessionAnimal = namesplit{2};
        disp(['Session loading: ',sessionName]);

        % Store session information in struct
        targetStruct(sessionID).mouse = sessionAnimal;
        targetStruct(sessionID).session = sessionName;
        targetStruct(sessionID).type = task;

        % Load recording data
        % Some session might not have LJ or NI photometry, initialize related
        % variable to 0 to avoid error
        rollingGreenLP = 0; timePhotometry = 0; photometryNI = 0;
        trials = 0; airpuff_rounded = 0;
        load(strcat(addlist{i},'\','sync_',sessionName,'.mat'),...
            'trials','params','rightLick','blueLaser',...
            'timeNI','timePhotometry','rollingGreenLP','photometryNI',...
            'rightSolenoid','airpuff_rounded','airpuff');

        % For old sessions
        if ~isfield(params.sync,'timeNI'); params.sync.timeNI = timeNI; end
        if ~isfield(params.sync,'timePhotometry'); params.sync.timePhotometry = timePhotometry; end
        if ~isfield(params.sync,'behaviorFs'); params.sync.behaviorFs = 10000; end
        if airpuff_rounded == 0; airpuff_rounded = airpuff; end
        
        % Get session cutoff
        [trials,cutoff_sample] = getSessionCutoff(trials,task,...
                            max_misses_allowed=options.maxMissesAllowed,...
                            windowlength=options.performanceEvalLength);
        params.analysis.cutoff_sample = cutoff_sample;

        % Store session recordings in struct
        targetStruct(sessionID).params = params;
        targetStruct(sessionID).trials = trials;
        targetStruct(sessionID).licks = find(rightLick);
        targetStruct(sessionID).water = find(rightSolenoid);
        targetStruct(sessionID).airpuff = find(airpuff_rounded);
        targetStruct(sessionID).stim = trials{trials.isTone == 0 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","performing"]};
        targetStruct(sessionID).tone = trials{trials.isTone == 1 & trials.isStim == 0,["TrialNumber","CueTime","OutcomeTime","performing"]};
        targetStruct(sessionID).pair = trials{trials.isTone == 1 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","performing"]};
        targetStruct(sessionID).randomShutter = round(rand([500,1])*length(blueLaser));
        targetStruct(sessionID).photometryLJ = rollingGreenLP;
        targetStruct(sessionID).photometryNI = photometryNI;

        clear dirsplit namesplit params rollingGreenLP photometryNI ...
            trials rightLick airpuff_rounded timeNI timePhotometry ...
            rightSolenoid 
        disp(['Session loaded: ',sessionName]);
    end
    
end

end