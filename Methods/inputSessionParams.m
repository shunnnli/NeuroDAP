function [Answer,Canceled] = inputSessionParams(sessionList,options)

arguments
    sessionList cell
    options.paradigm double = 1
    options.pavlovian logical = false
    options.reactionTime double = 2
    options.minLicks double = 2

    options.optoTriggered logical = false
    options.optoInverted logical = true
    options.optoPulseFreq double = 50
    options.optoPulseDuration double = 5
    options.optoStimDuration double = 500
end

Prompt = {'Session',...
          'Paradigm','Pavlovian',...
          'ReactionTime','minLicks',...
          'OptoTriggered','OptoInverted',...
          'OptoPulseFreq','OptoPulseDuration','OptoStimDuration'};

Formats(1,1).enable = 'inactive'; % SessionName, acts as a label
Formats(2,1).style = 'popupmenu'; % Behavior paradigm
Formats(2,1).type = 'list';
Formats(2,1).items = {'random','reward pairing','punish pairing'};
Formats(3,1).type = 'check';      % Pavlovian or operant
Formats(4,1).type = 'edit';       % Reaction time
Formats(5,1).type = 'edit';       % minLicks
Formats(6,1).type = 'check';      % Opto pulse signaled triggered or full event
Formats(7,1).type = 'check';      % Opto pulse inverted or not
Formats(8,1).type = 'edit';       % Opto pulse frequency
Formats(9,1).type = 'edit';       % Opto pulse duration
Formats(10,1).type = 'edit';      % Opto stim duration

% Create initial answer
DefAns = cell(size(Prompt,1),length(sessionList));
for s = 1:length(sessionList)
    % session name
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    DefAns{1,s} = sessionName(1:end-3);
    
    % Behavior params
    DefAns{2,s} = options.paradigm; % random
    DefAns{3,s} = options.pavlovian; % pavlovian or operant
    DefAns{4,s} = num2str(options.reactionTime);
    DefAns{5,s} = num2str(options.minLicks);

    % Opto params
    DefAns{6,s} = options.optoTriggered;
    DefAns{7,s} = options.optoInverted;
    DefAns{8,s} = num2str(options.optoPulseFreq);
    DefAns{9,s} = num2str(options.optoPulseDuration);
    DefAns{10,s} = num2str(options.optoStimDuration);
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
