
minAPs=1;

maxZones=getZone([]);
zFigure=[];
pulseAmp=-100;

for zz=1:maxZones
	zFigure(zz)=figure;
	title(['Sag Zone ' num2str(zz)]);
	hold on
	apResults.(['Z' num2str(zz)]).sagAvgTrace=[];
	apResults.(['Z' num2str(zz)]).sagV=[];	
	apResults.(['Z' num2str(zz)]).sagML=[];	
end

for counter=1:length(csAllCells)
	newCell=csAllCells(counter);
	disp([newCell.newName])
	aps=find((newCell.pulseI==pulseAmp) & (newCell.traceQC) & ~isnan(newCell.pulseV)); 
	
	if ~isempty(aps) 
		disp([num2str(counter) '     YES' num2str(zoneList)]);
		aps=aps(1);
		
		zoneList=getZone(newCell);
		
		if isempty(find(zoneList==1, 1))
			zoneList=[1 zoneList];
		end
		
		for zz=zoneList
			newData=newCell.acq{aps}.data;
			if isempty(newCell.Folder)
				newData=newCell.acq{aps}.data;
			else
				newData=[nan(1,1000) newCell.acq{aps}.data(1:7500) nan(1,4000) newCell.acq{aps}.data(7501:end)];
			end
			
			if isempty(apResults.(['Z' num2str(zz)]).sagAvgTrace)
				apResults.(['Z' num2str(zz)]).sagAvgTrace=newData;
			else
				apResults.(['Z' num2str(zz)]).sagAvgTrace=...
					apResults.(['Z' num2str(zz)]).sagAvgTrace+newData;
			end
						
			apResults.(['Z' num2str(zz)]).sagV(end+1)=newCell.sagV(aps);
%			apResults.(['Z' num2str(zz)]).sagML(end+1)=newCell.ML;

			figure(zFigure(zz));
			plot(newData);
		end
	end
end

for zz=1:maxZones
	apResults.(['Z' num2str(zz)]).sagAvgTrace=...
		apResults.(['Z' num2str(zz)]).sagAvgTrace / ...
		length(apResults.(['Z' num2str(zz)]).sagV);
	figure(zFigure(zz));
	[apResults.(['Z' num2str(zz)]).sagAvgTraceRecalc, apResults.(['Z' num2str(zz)]).sagAvgTraceSD]= ...
		phSummary_AvgSTDWindow(1);
	title(['Sag Zone ' num2str(zz) ' AVG +/- SD']);
end

