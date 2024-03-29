% Caluclate FDR significance for a given cluster, given two separate spike
% rates with the same sampling time
% event time: in samples
% Inputs: "spikeRate_for_each_cluster" from function "getSpikesByCluster": PxNxM
% matrix, spike rate for M neurons across N time points and P events; the cluster C you want to analyze
% Outputs: the matrix called "significance" which will be plugged in the other function "drawSpikeRateWithSignificance": a CxN significance matrix (0/1 matrix, indicating whether the FDR value is
% <0.05). Also give out a PxN FDR matrix (p value adjusted by multiple comparison) (MxN will take too long, so will just analyze interested clusters)

% ======================================= %
function [significance, pvalue] = ranksumSignificance(spikeRate1,spikeRate2,clusterList,options)

arguments
    spikeRate1 double
    spikeRate2 double
    clusterList double
    options.byTimeBin logical = false
end

% Q = ones(1,size(spikeRate1,2));
FDR = ones(length(clusterList),size(spikeRate1,3));
for i = 1:length(clusterList)
    neuron = clusterList(i);
    
    if options.byTimeBin
        for timebin = 1:size(spikeRate1,3)
            condition1 = spikeRate1(neuron, :, timebin);
            condition2 = spikeRate2(neuron, :, timebin);
            pvalue(i,timebin) = ranksum(condition1, condition2);
        end
%     %     [FDR, Q] = mafdr(pvalue(i,:)'); % mafdr has to contain all p-values for all clusters, pvalue has to be a column vector
%         disp(['Calculated significance for Cluster ',num2str(neuron)]);
    else
        condition1 = squeeze(mean(spikeRate1(neuron, :, :), 3));
        condition2 = squeeze(mean(spikeRate2(neuron, :, :), 3));
        if mean(condition1)==0 && mean(condition2)==0
            pvalue(i) = 1; 
        else
            pvalue(i) = ranksum(condition1, condition2);
        end
    end

end
% significance = (Q < 0.025);
significance = (pvalue < (0.05/size(spikeRate1,3)));  % bonferroni correction
% significance = (pvalue < 0.025);

end
