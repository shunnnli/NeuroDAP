allFields=fieldnames(apResults.Z1);
maxZones=getZone([]);
allSummary=[];

for fc=1:length(allFields)
	fns=allFields{fc};
	
	if length(apResults.Z1.(fns))<1000 % assume that long things are average traces
												% is you have more than
												% 1000 cells change this
												% number!
		disp(fns)

		ss=zeros(maxZones,3);
		for cc=1:maxZones
			ss(cc,:)=[length(apResults.(['Z' num2str(cc)]).(fns)) mean(apResults.(['Z' num2str(cc)]).(fns)) std(apResults.(['Z' num2str(cc)]).(fns))];
			disp(ss(cc,:))
		end
		disp('')

		allSummary.(fns)=ss;
		figure
		hold on
		bar(1:maxZones,ss(:,1))
		errorbar(1:maxZones,ss(:,1),ss(:,2),'.')
		title(fns)

		
	end
end

