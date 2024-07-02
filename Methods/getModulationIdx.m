function index = getModulationIdx(pre,post,clusterList,options)

arguments
    pre double
    post double
    clusterList double
    options.byTrials logical = false
end

index = nan(size(clusterList));
for i = 1:length(clusterList)
    neuron = clusterList(i);
    if options.byTrials
    
    else
        f_post = mean(squeeze(mean(post(neuron,:,:),2)));
        f_pre = mean(squeeze(mean(pre(neuron,:,:),2)));

        index(i) = (f_post - f_pre) / (f_post + f_pre);
    end
end

end