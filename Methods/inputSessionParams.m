function [Answer,Canceled] = inputSessionParams(sessionList,options)

arguments
    sessionList cell
    options.paradigm double = 1
    options.redStim logical = true
    options.pavlovian logical = false
    options.reactionTime double = 2
    options.minLicks double = 2

    options.optoTriggered logical = false
    options.optoInverted logical = true
    options.redPulseFreq double = 50
    options.redPulseDuration double = 5
    options.redStimDuration double = 500
    options.bluePulseFreq double = 30
    options.bluePulseDuration double = 10
    options.blueStimDuration double = 500

    options.includeOtherStim logical = true
end

Prompt = {'Session',...
          'Paradigm','redStim','Pavlovian',...
          'ReactionTime','minLicks',...
          'OptoTriggered','OptoInverted',...
          'RedPulseFreq','RedPulseDuration','RedStimDuration',...
          'BluePulseFreq','BluePulseDuration','BlueStimDuration',...
          'IncludeOtherStim'};

Formats(1,1).enable = 'inactive'; % SessionName, acts as a label
Formats(2,1).style = 'popupmenu'; % Behavior paradigm
Formats(2,1).type = 'list';
Formats(2,1).items = {'random','reward pairing','punish pairing'};
Formats(3,1).type = 'check';      % red stim or blue
Formats(4,1).type = 'check';      % Pavlovian or operant
Formats(5,1).type = 'edit';       % Reaction time
Formats(6,1).type = 'edit';       % minLicks
Formats(7,1).type = 'check';      % Opto pulse signaled triggered or full event
Formats(8,1).type = 'check';      % Opto pulse inverted or not
Formats(9,1).type = 'edit';       % Red pulse frequency
Formats(10,1).type = 'edit';       % Red pulse duration
Formats(11,1).type = 'edit';      % Red stim duration
Formats(12,1).type = 'edit';       % Blue pulse frequency
Formats(13,1).type = 'edit';       % Blue pulse duration
Formats(14,1).type = 'edit';      % Blue stim duration
Formats(15,1).type = 'check';       % Include other stim

% Create initial answer
DefAns = cell(size(Prompt,1),length(sessionList));
for s = 1:length(sessionList)
    % session name
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    DefAns{1,s} = sessionName(1:end-3);
    
    % Behavior params
    DefAns{2,s} = options.paradigm;  % random
    DefAns{3,s} = options.redStim;   % redStim
    DefAns{4,s} = options.pavlovian; % pavlovian or operant
    DefAns{5,s} = num2str(options.reactionTime);
    DefAns{6,s} = num2str(options.minLicks);

    % Opto params
    DefAns{7,s} = options.optoTriggered;
    DefAns{8,s} = options.optoInverted;
    DefAns{9,s} = num2str(options.redPulseFreq);
    DefAns{10,s} = num2str(options.redPulseDuration);
    DefAns{11,s} = num2str(options.redStimDuration);
    DefAns{12,s} = num2str(options.bluePulseFreq);
    DefAns{13,s} = num2str(options.bluePulseDuration);
    DefAns{14,s} = num2str(options.blueStimDuration);
    DefAns{15,s} = options.includeOtherStim;
end
DefAns = cell2struct(DefAns,Prompt,1);
Prompt = repmat(Prompt',1,2);

% Create inputsdlg
Title = 'Enter session and analysis params';
Options.AlignControls = 'on';
Options.CreateFcn = @(~,~,handles)celldisp(get(handles,'type'));
Options.DeleteFcn = @(~,~,handles)celldisp(get(handles,'type'));
[Answer,Canceled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

end
