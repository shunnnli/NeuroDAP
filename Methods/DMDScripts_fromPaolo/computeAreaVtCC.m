function [areaVtPulse, areaVtControl] = computeAreaVtCC(DataStruct, optoParameters)

   data = DataStruct.data;
   data = preprocessSignalVC(data);
   
   baselineStart = 5000;
   baselineEnd = min([15000, optoParameters.delayCurrentPulse, optoParameters.delayLightPulse]);
   baselineAverage = mean(data(baselineStart:baselineEnd));
   
   data = data - baselineAverage;
   
   widthCurrentPulse = optoParameters.widthCurrentPulse;
   currentPulseStart = optoParameters.delayCurrentPulse;
   currentPulseEnd = currentPulseStart + widthCurrentPulse;
   
   if currentPulseStart == 0
       
       widthCurrentPulse = 6000;
       currentPulseStart = 22000;
       currentPulseEnd = currentPulseStart + widthCurrentPulse;
       
   end
   
   areaVtPulse = trapz(currentPulseStart:currentPulseEnd,data(currentPulseStart:currentPulseEnd)) / 10000;

   areaVtControlEnd = size(data,2);
   areaVtControlStart = areaVtControlEnd - widthCurrentPulse;
   areaVtControl = trapz(areaVtControlStart:areaVtControlEnd,data(areaVtControlStart:areaVtControlEnd)) / 10000;

end

