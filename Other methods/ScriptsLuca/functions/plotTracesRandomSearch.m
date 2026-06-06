function plotTracesRandomSearch(wName, responseSign)

    global state

    % response and control window parameters for plotting
    durationPulseWindow = 500; % 50 ms
    durationControlWindow = durationPulseWindow;

    saveFolder = state.files.savePath;

    waveData = getWave(wName, 'data');
    hString = getfield(get(wName, 'UserData'), 'headerString');
    aiChanStr = state.phys.internal.lastLinesUsed;
    pulseString = valueFromHeaderString(['state.phys.internal.pulseString_' aiChanStr{2}], hString);
    
    delayFirstPulse = phUtil_parsePulsePatternString(pulseString, 'delay');
    isi = phUtil_parsePulsePatternString(pulseString, 'isi');
    nPulses = phUtil_parsePulsePatternString(pulseString, 'numPulses');
    pulsesTimeArray = (delayFirstPulse + linspace(0,(nPulses - 1)*isi,nPulses))*10;

    pulseTraces = cell(1, length(pulsesTimeArray));
    controlTraces = cell(1, length(pulsesTimeArray));
    pulsePeakData = nan(length(pulsesTimeArray),2);
    controlPeakData = nan(length(pulsesTimeArray),2);
    
    % No overlap in case standard RC parameters are used
    timeAfterEvent = 5000;
    startFirstBaseline = min(timeAfterEvent,pulsesTimeArray(1) - 1000);
    endFirstBaseline = pulsesTimeArray(1);
    startSecondBaseline = min(pulsesTimeArray(end) + timeAfterEvent,length(waveData) - 1000);
    endSecondBaseline = length(waveData);
    
    % Average and SD of combined baseline windows
    baselineAverage = mean([waveData(startFirstBaseline:endFirstBaseline),waveData(startSecondBaseline:endSecondBaseline)]);
    baselineSD = std([waveData(startFirstBaseline:endFirstBaseline),waveData(startSecondBaseline:endSecondBaseline)]);
    responseThreshold = 5*baselineSD;
    
    figure;

    for iPulse = 1:nPulses

        pulseStart = round(pulsesTimeArray(iPulse));
        pulseEnd = pulseStart + durationPulseWindow;
        
        controlStart = pulseStart - durationControlWindow;
        controlEnd = pulseStart;
        
        pulseTraces{iPulse} = waveData(pulseStart:pulseEnd) - baselineAverage;
        pulseTraces{iPulse} = preprocessSignalVC(pulseTraces{iPulse});
        
        controlTraces{iPulse} = waveData(controlStart:controlEnd) - baselineAverage;
        controlTraces{iPulse} = preprocessSignalVC(controlTraces{iPulse});
        
        pulsePeakData(iPulse,:) = extractPeakData(pulseTraces{iPulse}, responseThreshold, responseSign);
        controlPeakData(iPulse,:) = extractPeakData(controlTraces{iPulse}, responseThreshold, responseSign);
        
        plot(pulseTraces{iPulse},'Color',[0 0.4470 0.7410],'LineWidth',0.1);
        plot(controlTraces{iPulse},'Color',[0.4 0.4 0.4],'LineWidth',0.1);
        
    end
    
    hline = yline(responseSign*responseThreshold+baselineAverage, 'r--', 'LineWidth', 1);
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    meanPulseTraces = mean(cell2mat(pulseTraces(:)'), 2);
    meanControlTraces = mean(cell2mat(controlTraces(:)'), 2);
    
    plot(meanPulseTraces,'Color',[0 0.4470 0.7410],'LineWidth',0.6);
    plot(meanControlTraces,'Color',[0.4 0.4 0.4],'LineWidth',0.6);
    
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
    legend('Pulse','Control','Location', 'northeast', 'FontSize',12);
    
    figure;
    nBins = 10;
    
    h = histogram(pulsePeakData(:,1),nBins);
    h.FaceColor = [0 0.4470 0.7410];
    hold on
    
    h = histogram(controlPeakData(:,1),nBins);
    h.FaceColor = [0.4 0.4 0.4];
    
    set(gca, 'LineWidth', 0.5, 'FontSize', 12, 'FontName', 'Arial');
    legend('Pulse','Control','Location', 'northeast', 'FontSize',12);

end