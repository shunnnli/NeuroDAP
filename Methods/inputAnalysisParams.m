function [Answer,Canceled] = inputAnalysisParams(sessionList,options)

arguments
    sessionList cell
    options.reloadAll logical = false
    options.recordLJ  char = '[1,1,0]'
    options.rollingWindowTime double = 180
    options.withPhotometryNI logical = false
    options.plotPhotometry logical = true
    options.plotBehavior logical = true
end

Prompt = {'Session','reloadAll',...
          'recordLJ','rollingWindowTime',...
          'withPhotometryNI',...
          'plotPhotometry','plotBehavior'};

Formats(1,1).enable = 'inactive'; % SessionName, acts as a label
Formats(2,1).type = 'check';
Formats(3,1).type = 'edit';       % labjack.record
Formats(4,1).type = 'edit';       % rollingWindowTime
Formats(5,1).type = 'check';      % With photometryNI
Formats(6,1).type = 'check';      % Plot photometry
Formats(7,1).type = 'check';      % Plot behavior

% Create initial answer
DefAns = cell(size(Prompt,1),length(sessionList));
for s = 1:length(sessionList)
    % session name
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    DefAns{1,s} = sessionName(1:end-3);
    
    % Behavior params
    DefAns{2,s} = options.reloadAll;
    DefAns{3,s} = options.recordLJ;
    DefAns{4,s} = num2str(options.rollingWindowTime);
    DefAns{5,s} = options.withPhotometryNI;
    DefAns{6,s} = options.plotPhotometry;
    DefAns{7,s} = options.plotBehavior;
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