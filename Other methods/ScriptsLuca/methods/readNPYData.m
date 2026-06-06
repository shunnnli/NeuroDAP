function varargout = readNPYData(ksDir)

% Read spike_times.npy, spike_clusters.npy
% Return    clustergroup = [cluster_id group KSLabel]
%           spike_times
%           spike_clusters
% Need to install npy-matlab

% Load cluster_KSLabel.tsv
clusterpath = fullfile(ksDir,'cluster_KSLabel.tsv');
fid = fopen(clusterpath);
C = textscan(fid,'%d%C','HeaderLines',1);
fclose(fid);
cluster_id = C{1,1};
KSLabel = C{1,2};
ksgroup = table(cluster_id,KSLabel);

% Load session cluster_group.tsv
clusterpath = fullfile(ksDir,'cluster_group.tsv');
fid = fopen(clusterpath);
C = textscan(fid,'%d%C','HeaderLines',1);
fclose(fid);

cluster_id = C{1,1}; % cluster_id that I labeled on phy
group = C{1,2};
sorted = table(cluster_id,group);
cluster_id = setdiff(sorted.cluster_id,ksgroup.cluster_id);

% Add merged units to clustergroup (units without KSLabel)
KSLabel = repelem(categorical({'merged'}),length(cluster_id),1);
if ~isempty(cluster_id) && ~isempty(KSLabel)
    noKSLabel = table(cluster_id,KSLabel);
else
    noKSLabel = [];
end
ksgroup = [ksgroup;noKSLabel];

% Create label table label by phy
cluster_id = setdiff(ksgroup.cluster_id,sorted.cluster_id);
group = repelem(categorical({'unsorted'}),length(cluster_id),1);
unsorted = table(cluster_id,group);
labelgroup = [sorted;unsorted];

% Merge to form final table
clustergroup = join(labelgroup,ksgroup);
varargout{1} = clustergroup;

% Load spike_times, spike_clusters
spike_times = readNPY(strcat(ksDir,'/spike_times.npy'));
spike_clusters = readNPY(strcat(ksDir,'/spike_clusters.npy'));
varargout{2} = spike_times; varargout{3} = spike_clusters;

% Load cluster_info.tsv
clusterpath = fullfile(ksDir,'cluster_info.tsv');
try
    cluster_info = tsvread(clusterpath);
    cluster_info(1,:) = [];
    varargout{4} = cluster_info;
catch
    warning('cluster_info.tsv not found. Might need to run phy.');
    cluster_info = [];
    varargout{4} = cluster_info;
end


end