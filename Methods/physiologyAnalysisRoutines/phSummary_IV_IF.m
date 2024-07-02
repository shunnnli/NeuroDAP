%% set parameters

% use all amplitudes or only the ones listed below
useAllAmplitudes=1;
useFirstOnly=1; % if above is =1, use all the traces or only those on the first pass

% use the bins or the exact?
useExact=1;

% the list of pulse amplitudes to keep and analyze.  It doesn't matter if
% there are entries here that are not used
bExact=[-100 -50 0 20 40 50 60 80 100 150 200 250];
bExact=-100:50:900;

% instead define bins by low and high
bLow= [-100 -60 0  1 40  80 100 101 130 180 ];
bHigh=[-100 -40 0 39 79 99 100 119 179 219 ];

% what is the minimum number of APs the cell has to have fired across all
% trials to be included.  At least 1?
minAPs=1;

if useAllAmplitudes
	pCount=0;
	bExact=[];
else
	if useExact
		pCount=length(bExact);
	else
		pCount=length(bLow);
	end
end

%% initialize
maxZones=getZone([]); % edit the function getZone to make it return what zone a cell belongs to

avgV=zeros(maxZones,pCount);
avgI=zeros(maxZones,pCount);
avgA=zeros(maxZones,pCount);
avgN=zeros(maxZones,pCount);

stdV=zeros(maxZones,pCount);
stdI=zeros(maxZones,pCount);
stdA=zeros(maxZones,pCount);
stdN=zeros(maxZones,pCount);

pulseListForThisCell=[];

			
%% run through the cells
for runThru=1:2 % first run calcualtes the averages.  Second does the Standard Deviation
	if runThru==2 && useAllAmplitudes
		bExact=[];
	end
	
	for counter=1:length(csAllCells)
		newCell=csAllCells(counter);

		if (newCell.QC==1) && (sum(newCell.nAP, 'omitnan')>=minAPs) % the cell should be included 
			disp([newCell.newName])
			zone=getZone(newCell);
			if zone==1
				zoneList=1;
			else
				zoneList=[1 zone];
			end

			pulseListForThisCell=[];
			if useAllAmplitudes
				for traceCounter=find((newCell.traceQC==1) & ~isnan(newCell.pulseV))
					pulseI=newCell.pulseI(traceCounter);
					bCounter=find(bExact==pulseI);
					if isempty(bCounter) %never seen this pulse before
						bExact(end+1)=pulseI;
						bCounter=length(bExact);
						pulseListForThisCell(end+1)=pulseI;
					elseif useFirstOnly && any(pulseListForThisCell==pulseI)
						bCounter=[]; % we did it before - skip it
					end
					if ~isempty(bCounter) % we will analyze it
						pulseListForThisCell(bCounter)=pulseI;

						for ind=zoneList  % run through it once to put in the ALL pile and then in zone specific pile
							if runThru==1 % first time through calculate the average
								if bCounter>size(avgN,2)
									avgV(:, bCounter)=0;
									avgI(:, bCounter)=0;
									avgA(:, bCounter)=0;
									avgN(:, bCounter)=0;
								end
								avgV(ind, bCounter)=avgV(ind, bCounter)+newCell.pulseV(traceCounter);
								avgI(ind, bCounter)=avgI(ind, bCounter)+newCell.pulseI(traceCounter);
								avgA(ind, bCounter)=avgA(ind, bCounter)+newCell.nAP(traceCounter);
								avgN(ind, bCounter)=avgN(ind, bCounter)+1;
							else % first time through calculate the standard deviation
								if bCounter>size(stdN,2)
									stdV(:, bCounter)=0;
									stdI(:, bCounter)=0;
									stdA(:, bCounter)=0;
									stdN(:, bCounter)=0;
								end
								stdV(ind, bCounter)=stdV(ind, bCounter)+(newCell.pulseV(traceCounter)-avgV(ind, bCounter))^2;
								stdI(ind, bCounter)=stdI(ind, bCounter)+(newCell.pulseI(traceCounter)-avgI(ind, bCounter))^2;
								stdA(ind, bCounter)=stdA(ind, bCounter)+(newCell.nAP(traceCounter)-avgA(ind, bCounter))^2;
								stdN(ind, bCounter)=stdN(ind, bCounter)+1;
							end				
						end
					end
				end
			else
				for bCounter=1:pCount
					ff1=[];
					if useExact
						ff=find((newCell.pulseI==bExact(bCounter)));
						if ~isempty(ff)
							ff1=ff(1);
						end
					else
						ff=find((newCell.pulseI>=bLow(bCounter)) & (newCell.pulseI<=bHigh(bCounter)));
						if ~isempty(ff)
							ff1=ff(1);
						end
					end
					
					if ~isempty(ff1)  && (newCell.traceQC(ff1)==1) && ~isnan(newCell.pulseV(ff1)) % trace passes QC
						zone=getZone(newCell);

						if zone==1
							zoneList=1;
						else
							zoneList=[1 zone];
						end
					
						for ind=zoneList  % run through it once to put in the ALL pile and then in zone specific pile
							if runThru==1 % first time through calculate the average
								avgV(ind, bCounter)=avgV(ind, bCounter)+newCell.pulseV(ff1);
								avgI(ind, bCounter)=avgI(ind, bCounter)+newCell.pulseI(ff1);
								avgA(ind, bCounter)=avgA(ind, bCounter)+newCell.nAP(ff1);
								avgN(ind, bCounter)=avgN(ind, bCounter)+1;
							else % first time through calculate the standard deviation
								stdV(ind, bCounter)=stdV(ind, bCounter)+(newCell.pulseV(ff1)-avgV(ind, bCounter))^2;
								stdI(ind, bCounter)=stdI(ind, bCounter)+(newCell.pulseI(ff1)-avgI(ind, bCounter))^2;
								stdA(ind, bCounter)=stdA(ind, bCounter)+(newCell.nAP(ff1)-avgA(ind, bCounter))^2;
								stdN(ind, bCounter)=stdN(ind, bCounter)+1;
							end
						end
					end
				end
			end
		end
	end
	
	if runThru==1
		avgI=avgI./avgN;
		avgV=avgV./avgN;
		avgA=avgA./avgN;
	else
		stdI=(stdI./stdN).^(0.5);
		stdV=(stdV./stdN).^(0.5);
		stdA=(stdA./stdN).^(0.5);
		semI=stdI./(stdN.^(0.5));
		semV=stdV./(stdN.^(0.5));
		semA=stdA./(stdN.^(0.5));
	end
end

if useAllAmplitudes % resort the data
	[aa, ia]=sort(avgI(1,:));
	avgI=avgI(:, ia);
	avgV=avgV(:, ia);
	avgA=avgA(:, ia);
	avgN=avgN(:, ia);
	stdI=stdI(:, ia);
	stdV=stdV(:, ia);
	stdA=stdA(:, ia);
	stdN=stdN(:, ia);		
end

rr=1:23;

figure;
for c=1:2
	vv=avgV(c,:);
	ii=avgI(c,:);
	vv(isnan(vv))=[];
	ii(isnan(ii))=[];
	plot(ii(rr), vv(rr), 'k')
	hold on
	vvs=semV(c,:);
	iis=semI(c,:);
	vvs(isnan(vvs))=[];
	iis(isnan(iis))=[];
	hh=errorbar(ii(rr), vv(rr), vvs(rr), vvs(rr));
	set(hh, 'LineWidth', c, 'color', 'k')
	title(['I V']);
end

figure;
for c=1:1
	aa=avgA(c,:);
	ii=avgI(c,:);
	aa(isnan(aa))=[];
	ii(isnan(ii))=[];
	
	plot(ii(rr), aa(rr), 'color', 'k')
	hold on
	aas=semA(c,:);
	iis=semI(c,:);
	aas(isnan(aas))=[];
	iis(isnan(iis))=[];
	hh=errorbar(ii(rr), aa(rr), aas(rr), aas(rr));
	set(hh, 'LineWidth', c, 'color', 'k')
	title(['I F']);	
end

figure;
for c=1:maxZones
	aa=avgA(c,:);
	vv=avgV(c,:);
	aa(isnan(aa))=[];
	vv(isnan(vv))=[];
	
	plot(vv, aa)
	hold on
	aas=semA(c,:);
	vvs=semV(c,:);
	aas(isnan(aas))=[];
	vvs(isnan(vvs))=[];
	hh=errorbar(vv, aa, aas, aas);
	set(hh, 'LineWidth', c)
	title(['V F']);
end