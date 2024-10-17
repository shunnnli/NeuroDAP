function [Answer,Canceled] = inputSessionParams_singleSlice(sessionList,options)

arguments
    sessionList cell
    options.paradigm double = 1
    options.redStim logical = true
    options.animal char = 'SL'

    options.reload logical = false
    options.calculateQC logical = false
    options.timeRange char = '[-20,100]'
    options.nArtifactSamples char = '10'
end

Prompt = {'Session','Animal','reload','calculateQC','redStim','Paradigm','timeRange','nArtifactSamples'};

Formats(1,1).enable = 'inactive'; % SessionName, acts as a label
Formats(2,1).type = 'edit';       % Mouse ID
Formats(3,1).type = 'check';       % reload
Formats(4,1).type = 'check';      % calculateQC
Formats(5,1).type = 'check';      % red stim or blue during behavior
Formats(6,1).style = 'popupmenu'; % Behavior paradigm
Formats(6,1).type = 'list';
Formats(6,1).items = {'random','reward pairing','punish pairing','reward ctrl','punish ctrl'};
Formats(7,1).type = 'edit';       % timeRange
Formats(8,1).type = 'edit';       % nArtifactSamples


% Create initial answer
DefAns = cell(size(Prompt,1),length(sessionList));
for s = 1:length(sessionList)
    % session name
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    DefAns{1,s} = sessionName;
    DefAns{2,s} = options.animal; % mouseID
    DefAns{3,s} = options.reload;
    DefAns{4,s} = options.calculateQC;
    DefAns{5,s} = options.redStim;   % redStim
    DefAns{6,s} = options.paradigm;  % random
    DefAns{7,s} = options.timeRange;
    DefAns{8,s} = options.nArtifactSamples;
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
