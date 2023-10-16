% Same as function drawSpikeRate, but label significant points with different colors 
% Inputs:  NOT same as drawSpikeRate (spikeRate: MxN matrix for averaged firing rate, containing M
% neurons and N time bins; clusterList: the list of all neurons you want to analyze) +
% the significance value for a certain cluster (1xN)
% Outputs: plots with spike rate and significance

% USE MOVING AVERAGE OF 5 TIME BINS = 50 ms 
% Calculate standard error of mean (SEM)
% ======================================= %
% temp = single_cell_activity_pinch{1};
% temp_1 = nanmean(temp,2);
% temp_2 = nanstd(temp,0,2)./sqrt(sum(~isnan(temp(1,:))));
% h = area(-5:0.1:10,[(temp_1-temp_2),(2*temp_2)],-2);
% set(h(1),'Edgecolor','none','FaceColor','none');
% set(h(2),'Edgecolor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
% plot(-5:0.1:10,mean(single_cell_activity_pinch{1},2),'k');
% ylim([-1.5 2.5])

function [] = drawSpikeRateWithSignificance(eventName,timeRange,...
                spikeRateByCluster,clusterList,textOn,ap,colors,significance)

% Set colormap
cmap = colormap(flipud(jet(length(clusterList))));

% Plot individual traces
for i = 1:length(clusterList)
    neuron = clusterList(i);
    temp = spikeRateByCluster(:,:,neuron);
    temp_1 = nanmean(temp,1); % firing rate vector, same length as time

    if sum(temp_1) > 10
%         smoothed_fr = smoothdata(spikeRate(neuron,:),'movmean',5);
%         smoothed_fr = spikeRate(neuron,:);
        smoothed_fr = smoothdata(temp_1,'movmean',5);
        smoothed_fr_significant_only = smoothed_fr.*significance;
        smoothed_fr_significant_only(find(smoothed_fr_significant_only==0)) = NaN;

        [~,x_max] = max(abs(smoothed_fr));
        t = linspace(timeRange(1),timeRange(2),size(spikeRateByCluster,2));
        temp_2 = nanstd(temp,0,1)./sqrt(size(temp,1));%sum(~isnan(temp(2,:)))); SEM = SD/sqrt(N), N is observation

        if length(clusterList) > length(colors)
%             plot(t,smoothed_fr,'-','LineWidth',1); 
            h = area(t,[(smoothed_fr-temp_2);(2*temp_2)]',0);
            set(h(1),'Edgecolor','none','FaceColor','none');
            set(h(2),'Edgecolor','none','FaceColor',colors(i),'FaceAlpha',0.3);
            plot(t,smoothed_fr,'-','LineWidth',1);
            hold on;
            scatter(t, smoothed_fr_significant_only, 17,'filled'); hold on
        else
            h = area(t,[(smoothed_fr-temp_2);(2*temp_2)]',0);
            set(h(1),'Edgecolor','none','FaceColor','none');
            set(h(2),'Edgecolor','none','FaceColor',colors(i),'FaceAlpha',0.3);  
            plot(t,smoothed_fr,'-','LineWidth',1,'Color',colors(i));
            hold on;            
            scatter(t, smoothed_fr_significant_only, 17,'filled'); hold on
        end
        drawnow
        
        if textOn
            xtxtpos = t(x_max);
            ytxtpos = smoothed_fr(x_max);
            text(xtxtpos,ytxtpos,num2str(ap.goodClusters(neuron))); hold on
        end
    end
end

% Set plot related params
% set(gca,'ColorOrder',cmap);

ylim([0 inf])
xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,100]);
xline(0,':',eventName,'Color', [.6 .6 .6],'LineWidth',1.5);
box off

end % drawPSTH