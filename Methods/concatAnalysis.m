 function summary = concatAnalysis(sessionList)

arguments
    sessionList cell
end

summary = struct([]);
disp('Ongoing: summary.mat not found, creating a new one');

% Group sessions to summary
for s = 1:length(sessionList)
    % Find sessionName
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit

    % Load session
    load(strcat(sessionList{s},filesep,'analysis_',sessionName,'.mat'),'analysis');

    summary = [summary,analysis];
    disp(['Finished: session ', sessionName,' loaded (',...
          num2str(s),'/',num2str(length(sessionList)),')']);
end

% Check summary format (should all be chars NOT strings)
stringColumnsLabels = {'animal','date','session','task','event','name','system'};
for i = 1:length(stringColumnsLabels)
    for row = 1:length(summary)
        if isstring(summary(row).(stringColumnsLabels{i})) 
            summary(row).(stringColumnsLabels{i}) = convertStringsToChars(summary(row).(stringColumnsLabels{i}));
        end
    end
end

end