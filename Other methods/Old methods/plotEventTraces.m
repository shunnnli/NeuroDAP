function fig = plotEventTraces(eventIdx,baselineIdx,longTimeRange,shortTimeRange,groupSize,color)

binSize = params.finalTimeStep;

for i = 1:length(eventIdx)
    [~, eventTime_lj] = min(abs(timePhotometry-timeNI(eventIdx(i))));
    [~, baselineTime_lj] = min(abs(timePhotometry-timeNI(baselineIdx(i))));
    eventInLJ(i) = eventTime_lj;
    baselineInLJ(i) = baselineTime_lj;
end

fig = initializeFig(0.5,1);
subplot(4,1,1)
[traces,t] = getTraces(eventInLJ/params.finalFs,rollingGreenLP,longTimeRange,binSize);
[baseline,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,longTimeRange,binSize);
plotCI(t,baseline(1:end,:),[.75 .75 .75]);
plotCI(t,traces(1:end,:),color(1,:));
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label); 
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    [label,' (n=',num2str(length(eventIdx)),')']},...
    'Location','northeast'); 

subplot(4,1,3)
nLines = ceil(size(traces,1)/groupSize);
legendList = cell(nLines,1);
nColors = round(linspace(1,size(color,1),nLines));
for i = 1:nLines
    startTrial = (i-1)*groupSize+1; 
    if i == nLines; endTrial = size(traces,1);
    else; endTrial = i*groupSize; end
    plotCI(t,traces(startTrial:endTrial,:),color(nColors(i),:));
    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
end
plotEvent(label,eventDuration,'r');
legend(legendList);

subplot(4,1,2)
[traces,t] = getTraces(eventInLJ/params.finalFs,rollingGreenLP,shortTimeRange,binSize);
[baseline,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,shortTimeRange,binSize);
plotCI(t,baseline(1:end,:),[.75 .75 .75]);
plotCI(t,traces(1:end,:),color(1,:));
plotEvent(label,eventDuration,'r'); 
xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    [label,' (n=',num2str(length(eventIdx)),')']},...
    'Location','northeast'); 

subplot(4,1,4)
nLines = ceil(size(traces,1)/groupSize);
legendList = cell(nLines,1);
nColors = round(linspace(1,size(color,1),nLines));
for i = 1:nLines
    startTrial = (i-1)*groupSize+1; 
    if i == nLines; endTrial = size(traces,1);
    else; endTrial = i*groupSize; end
    plotCI(t,traces(startTrial:endTrial,:),color(nColors(i),:));
    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
end
plotEvent(label,eventDuration,'r');
legend(legendList);

saveas(gcf,strcat(session.path,'\psth_',label,'_',session.name,'.png'));

end