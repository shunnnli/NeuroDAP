function trialTables = loadTrialTables(sessionList,options)

arguments
    sessionList cell
    options.animalIdxInName double = 2 % where animal name will appear in sessionName
end

data = {}; name = {}; 
animalList = {};
sessionIdx = [];

for s = 1:length(sessionList)
    % Find animal name and session name
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; 

    dirsplit = strsplit(sessionName,{'-','_'});
    animalName = dirsplit{options.animalIdxInName};

    % When there is a new animal
    if ~any(contains(animalList,animalName))
        animalList = [animalList,animalName];
        sessionIdx(end+1) = 1;
        animalIdx = length(animalList);
    else
        animalIdx = find(contains(animalList,animalName));
        sessionIdx(animalIdx) = sessionIdx(animalIdx) + 1;
    end
    
    % Load trial table
    load(strcat(sessionList{s},filesep,'behavior_',sessionName,'.mat'),'trials');
    data{animalIdx,sessionIdx(animalIdx)} = trials;
    name{animalIdx,sessionIdx(animalIdx)} = sessionName;
    disp(['Finished: trial table for ',animalName,' session ',sessionName, ' loaded.']);
end

trialTables.data = data;
trialTables.name = name;

end