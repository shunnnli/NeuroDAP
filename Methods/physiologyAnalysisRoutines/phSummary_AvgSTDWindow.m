function [ meanData, stdData ] = phSummary_AvgSTDWindow(doPlot)
% Click on a window and run this.  It returns the avg and the std of all
% the data in the plots

	if nargin<1
		doPlot=0;
	end
	

	ll=get(gca, 'Children');
	dCount=length(ll);
	for counter=1:dCount
		newData=ll(counter).YData;
		if counter==1
			meanData=newData/dCount;
		else
			meanData=meanData+newData/dCount;
		end
	end
	stdData=meanData.*0;
	for counter=1:dCount
		newData=ll(counter).YData;
		stdData=stdData+(newData-meanData).^2/dCount;
	end
	stdData=sqrt(stdData)/sqrt(dCount);
	
	if doPlot
		figure
		hold on
		shadedErrorBar(1:length(meanData), meanData, stdData)
	end
end

