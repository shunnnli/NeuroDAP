%% Shun_EPLHbsimulation

% Initialize EP-LHb synapse matrix
nLHb = 10;
nSynapses = 100;
nSimulations = 1000;

totalEIsum_peaks = []; totalEIsum_aucs = []; % 10000 data points
totalEIindex_peaks = []; totalEIindex_aucs = []; % 10000
totalEP_peaks = []; totalEP_aucs = []; % 1000

for i = 1:nSimulations
    % Step 1: randomly select 10 cells from baseline
    LHbCellsIdx = randomIdx(randperm(length(randomIdx), 10));
    totalEPSC_aucs = EPSC_aucs(LHbCellsIdx);
    totalIPSC_aucs = IPSC_aucs(LHbCellsIdx);
    totalEPSC_peaks = EPSC_peaks(LHbCellsIdx);
    totalIPSC_peaks = IPSC_peaks(LHbCellsIdx);
    
    % Step 2: determine synapse strength by randomly select from a normal
    % distribution from 0 to 1
    weights = abs(randn(nLHb,nSynapses));
    weights = weights ./ sum(weights,2);
    synapseEPSC_peaks = weights .* totalEPSC_peaks;
    synapseIPSC_peaks = weights .* totalIPSC_peaks;
    syanpseEPSC_aucs = weights .* totalEPSC_aucs;
    synapseIPSC_aucs = weights .* totalIPSC_aucs;

    % Calculate EI index
    curEIindex_peaks = abs(totalIPSC_peaks)-abs(totalEPSC_peaks)./abs(totalIPSC_peaks)+abs(totalEPSC_peaks);
    curEIindex_aucs = abs(totalIPSC_aucs)-abs(totalEPSC_aucs)./abs(totalIPSC_aucs)+abs(totalEPSC_aucs);

    % Store statistics
    totalEIsum_peaks = [totalEIsum_peaks; sum(synapseEPSC_peaks,2) + sum(synapseIPSC_peaks,2)];
    totalEIsum_aucs = [totalEIsum_aucs; sum(syanpseEPSC_aucs,2) + sum(synapseIPSC_aucs,2)];
    totalEIindex_peaks = [totalEIindex_peaks; curEIindex_peaks];
    totalEIindex_aucs = [totalEIindex_aucs; curEIindex_aucs];
    totalEP_peaks = [totalEP_peaks; sum(synapseEPSC_peaks,'all') + sum(synapseIPSC_peaks,'all')];
    totalEP_aucs = [totalEP_aucs; sum(syanpseEPSC_aucs,'all') + sum(synapseIPSC_aucs,'all')];
end

%% Plot initialize distribution

close all;

initializeFig(1,1); master=tiledlayout(2,3);

nexttile(master,1);
plotDistribution(totalEP_peaks,tile=1,masterlayout=master,xlabel='Total EP amplitude across all LHb');

nexttile(master,4);
plotDistribution(totalEP_aucs,tile=4,masterlayout=master,xlabel='Total EP charge across all LHb');

nexttile(master,2);
plotDistribution(totalEIsum_peaks,tile=2,masterlayout=master,xlabel='EPSC+IPSC amplitude');

nexttile(master,5);
plotDistribution(totalEIsum_aucs,tile=5,masterlayout=master,xlabel='EPSC+IPSC charge');

nexttile(master,3);
plotDistribution(totalEIindex_peaks,tile=3,masterlayout=master,xlabel='EI amplitude index');

nexttile(master,6);
plotDistribution(totalEIindex_aucs,tile=6,masterlayout=master,xlabel='EI charge index');

%% Change weights

pctFlipped = 10;


%% code

initializeFig(1,1); tiledlayout(1,3);
nexttile;
imagesc(synapseEPSC_peaks); 
colormap("parula"); colorbar;
xlabel('Synapses'); ylabel('LHb neuron');

nexttile;
imagesc(synapseIPSC_peaks); colorbar;
xlabel('Synapses'); ylabel('LHb neuron');

nexttile;
imagesc(synapseEPSC_peaks + synapseIPSC_peaks); colorbar;
xlabel('Synapses'); ylabel('LHb neuron');