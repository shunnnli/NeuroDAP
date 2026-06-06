function [maxPSC, minPSC, meanPSC]=phAnalyzePSC(dataWave, stimAOchannel, options)
arguments
    dataWave struct
    stimAOchannel float
    options.vhold = -70
    options.vrev = 0
    options.analysisDelay = 5;
    options.analysisWindow = 20;
end

hString=dataWave.UserData.headerString;
deltaT=dataWave.xscale(2);

stimPattern=phUtil_HeaderValue(hString, ['state.phys.internal.pulseString_ao' num2str(stimAOchannel) ]);


stimNumber=phUtil_parsePulsePatternString(stimPattern, 'numPulses');
stimDelay=phUtil_parsePulsePatternString(stimPattern, 'delay');
stimISI=phUtil_parsePulsePatternString(stimPattern, 'isi');

startBaseline=max(1, (stimDelay-50)/deltaT);
endBaseline=(stimDelay-1)/deltaT;
baseline=mean(dataWave.data(startBaseline:endBaseline));

maxPSC=zeros(1, stimNumber);
minPSC=maxPSC;
meanPSC=maxPSC;

for pulseCounter=1:stimNumber

    startAna=stimDelay+(stimNumber-1)*stimISI+options.analysisDelay;
    endAna=startAna+options.analysisWindow;
    startAnaPt=1+startAna/deltaT;
    endAnaPt=1+endAna/deltaT;

    maxPSC(stimNumber)=max(dataWave.data(startAnaPt:endAnaPt))-baseline;
    minPSC(stimNumber)=min(dataWave.data(startAnaPt:endAnaPt))-baseline;
    meanPSC(stimNumber)=mean(dataWave.data(startAnaPt:endAnaPt)-baseline);
end



end






