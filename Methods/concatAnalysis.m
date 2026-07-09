function summary = concatAnalysis(sessionList, options)

arguments
    sessionList cell
    options.columnLabels = {'animal','date','session','task','event','name','system'}
    options.skipCamera logical = false
end

summary = struct([]);
errorSessions = {};
disp('Ongoing: summary.mat not found, creating a new one');

% Group sessions to summary
for s = 1:length(sessionList)
    % Find sessionName
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit

    % Load session
    load(strcat(sessionList{s},filesep,'analysis_',sessionName,'.mat'),'analysis');

    if exist('analysis','var')
        try
            % Skip camera rows if requested
            if options.skipCamera && ~isempty(analysis)
                systems = strings(1,length(analysis));
        
                for row = 1:length(analysis)
                    systems(row) = string(analysis(row).system);
                end
        
                analysis = analysis(systems ~= "Cam");
            end

            summary = [summary,analysis]; 

        catch e
            disp(strcat('  concatAnalysis: ', getReport(e)));
            errorSessions{end+1} = sessionList{s};
            continue
        end
    end
    disp(['Finished: session ', sessionName,' loaded (',...
          num2str(s),'/',num2str(length(sessionList)),')']);
end

% Check summary format (should all be chars NOT strings)
stringColumnsLabels = options.columnLabels;
for i = 1:length(stringColumnsLabels)
    for row = 1:length(summary)
        if isstring(summary(row).(stringColumnsLabels{i})) 
            summary(row).(stringColumnsLabels{i}) = convertStringsToChars(summary(row).(stringColumnsLabels{i}));
        end
    end
end

end