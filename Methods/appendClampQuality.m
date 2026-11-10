function combined = appendClampQuality(filename,current)
%APPENDCLAMPQUALITY Append new animals and replace previously saved animals.
% All tables for an updated animal are replaced together. Trials from the
% same animal in different datasets are not pooled by this function.
% Older files containing only clampQuality are backed up and upgraded.
% Missing trial/event details are left empty; unknown settings are recorded
% per animal in clampQualitySources.SettingsVerified.
% A changed responseWindow starts a fresh file after backing up all old
% results. Recalculate the other animals before appending them to this file.

arguments
    filename (1,1) string
    current (1,1) struct
end

required = {'clampQuality','clampQualityTrials', ...
            'clampQualityEventSummary','clampQualitySettings'};
if ~all(isfield(current,required))
    error('ClampQuality:GenerateFirst','Run section 1 (Generate clampQuality) first.');
end
if isempty(current.clampQuality)
    error('ClampQuality:NoData','No matching animals were generated; the saved file was not changed.');
end

settings = current.clampQualitySettings;
if ~isfield(settings,'sourceFiles'); settings.sourceFiles = strings(0,1); end
newAnimals = string(current.clampQuality.Animal);
sources = table(newAnimals,repmat({settings.sourceFiles},numel(newAnimals),1), ...
    repmat(settings.createdAt,numel(newAnimals),1),true(numel(newAnimals),1), ...
    'VariableNames',{'Animal','SourceFiles','GeneratedAt','SettingsVerified'});
combined = current;
legacyFields = {};
windowChanged = false;

if isfile(filename)
    previous = load(filename);
    if ~isfield(previous,'clampQuality') || ~istable(previous.clampQuality)
        error('ClampQuality:FileFormat', ...
            'The saved file must contain a clampQuality table; it was not changed.');
    end
    legacyFields = required(~isfield(previous,required));
    for k = 2:3
        field = required{k};
        if ~isfield(previous,field)
            % Ratios cannot recover the original trials or event statistics.
            emptyTable = current.(field);
            previous.(field) = emptyTable([],:);
        end
    end
    if ~isfield(previous,'clampQualitySettings')
        previous.clampQualitySettings = struct();
    end
    % Prevent appending measurements made with different definitions.
    compareFields = {'version','signal','events','baselineWindow', ...
        'responseWindow','performingOnly','requireCompleteWindows', ...
        'peakDefinition','pooling'};
    previousSettingsVerified = all(isfield(previous.clampQualitySettings,compareFields));
    for k = 1:numel(compareFields)
        field = compareFields{k};
        if ~isfield(previous.clampQualitySettings,field)
            continue; % unavailable in older exports; do not invent a value
        end
        if ~isfield(settings,field) || ...
                ~isequaln(previous.clampQualitySettings.(field),settings.(field))
            if strcmp(field,'responseWindow') && isfield(settings,field)
                windowChanged = true;
            else
                error('ClampQuality:SettingsMismatch', ...
                    'Analysis setting "%s" differs from the saved file; it was not changed.',field);
            end
        end
    end
    if ~windowChanged
        for k = 1:3
            field = required{k};
            oldTable = previous.(field);
            keep = ~ismember(lower(string(oldTable.Animal)),lower(newAnimals));
            combined.(field) = [oldTable(keep,:); current.(field)];
        end
        if isfield(previous,'clampQualitySources')
            oldSources = previous.clampQualitySources;
        else
            % Support the earlier single-export MAT format.
            oldAnimals = string(previous.clampQuality.Animal);
            oldFiles = strings(0,1);
            if isfield(previous.clampQualitySettings,'sourceFiles')
                oldFiles = previous.clampQualitySettings.sourceFiles;
            end
            oldCreatedAt = NaT;
            if isfield(previous.clampQualitySettings,'createdAt')
                oldCreatedAt = previous.clampQualitySettings.createdAt;
            end
            oldSources = table(oldAnimals,repmat({oldFiles},numel(oldAnimals),1), ...
                repmat(oldCreatedAt,numel(oldAnimals),1), ...
                repmat(previousSettingsVerified,numel(oldAnimals),1), ...
                'VariableNames',sources.Properties.VariableNames);
        end
        if ~ismember('SettingsVerified',oldSources.Properties.VariableNames)
            oldSources.SettingsVerified = repmat(previousSettingsVerified,height(oldSources),1);
        end
        keep = ~ismember(lower(string(oldSources.Animal)),lower(newAnimals));
        sources = [oldSources(keep,:); sources];
    end
end

combined.clampQuality = sortrows(combined.clampQuality,'Animal');
combined.clampQualityEventSummary = sortrows(combined.clampQualityEventSummary, ...
    {'Animal','Event','Condition'});
combined.clampQualityTrials = sortrows(combined.clampQualityTrials, ...
    {'Animal','Event','Condition','AnimalRow','TrialRow'});
combined.clampQualitySources = sortrows(sources,'Animal');
% Global settings describe verified data; retain an explicit list of animals
% imported without enough metadata to check their analysis definitions.
settings.unverifiedAnimals = string(sources.Animal(~sources.SettingsVerified));
settings.sourceFiles = strings(0,1);
for k = 1:height(sources)
    settings.sourceFiles = union(settings.sourceFiles,string(sources.SourceFiles{k}(:)));
end
combined.clampQualitySettings = settings;
[combined.clampQualityWeighted,combined.clampQualityWeights] = weightClampQuality( ...
    combined.clampQualityEventSummary,string(combined.clampQuality.Animal),string(settings.events));
combined.clampQualitySettings.weighting = ...
    'Common event weights = min(NClamp,NUnclamp); ratio of weighted means';

% Finish writing the new file before replacing any existing saved results.
folder = fileparts(filename);
if strlength(folder)==0; folder = string(pwd); end
if windowChanged || ~isempty(legacyFields)
    [~,baseName] = fileparts(filename);
    stamp = string(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    backupReason = "legacy";
    if windowChanged; backupReason = "previousWindow"; end
    backupFile = fullfile(folder,baseName + "_" + backupReason + "_" + stamp + ".mat");
    [ok,message] = copyfile(filename,backupFile);
    if ~ok
        error('ClampQuality:BackupFailed','Could not back up the older file: %s',message);
    end
end
temporaryFile = string(tempname(folder)) + ".mat";
save(temporaryFile,'-struct','combined');
[ok,message] = movefile(temporaryFile,filename,'f');
if ~ok
    delete(temporaryFile);
    error('ClampQuality:SaveFailed','Could not save clamp quality: %s',message);
end
if windowChanged
    warning('ClampQuality:ResponseWindowChanged', ...
        ['Response window changed from %s to %s seconds. Old results backed up to %s. ' ...
         'Started a new file with the %d currently generated animals. ' ...
         'Regenerate and append the other animals using the new window.'], ...
        mat2str(previous.clampQualitySettings.responseWindow), ...
        mat2str(settings.responseWindow),backupFile,numel(newAnimals));
elseif ~isempty(legacyFields)
    warning('ClampQuality:LegacyFile', ...
        ['Upgraded the older file; original backed up to %s. Missing fields: %s. ' ...
         'Trial/event details remain unavailable for animals not regenerated. ' ...
         'See clampQualitySettings.unverifiedAnimals for animals with unknown settings.'], ...
        backupFile,strjoin(legacyFields,', '));
end
end
