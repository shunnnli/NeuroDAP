% coreleaseAnalysis_Batch.m --- obsolete!
% v.1 by SAK 12.10.19 implemented batch mode of pre-existing code in varAroundSpots.m
clear all

% parameters:
cpaWind= 0; % in 0.1ms


% select main folder where we can find subfolder names & excel file
[~,mainpath] = uiputfile('*.*','Select main data folder', 'mainpath.mat');
cd(mainpath)
% Read excel sheet for getting acq #
LoadExcel('coreleaseExamples.xlsx')

for cell_i= 1:size(batchTableRaw,1)
 
clearvars -except batchTableRaw cpaWind cell_i mainpath
    
% pathlocation
FolderPath = strcat(mainpath,batchTableRaw{cell_i,2},'/Cell',string(batchTableRaw{cell_i,1}));
% move to folder
cd(FolderPath)
% define epoch number
epoch = batchTableRaw{cell_i,3}; 

% loads average trace and saves trial numbers
load(strcat(['AD0_e',num2str(epoch),'p2avg']))
traceNums = evalin('base',['AD0_e',num2str(epoch),'p2avg.UserData.Components']);

% saves pattern sequences
patternSeq=str2num(valueFromHeaderString('state.DMD.patternsString',...
  evalin('base',['AD0_e',num2str(epoch),'p2avg.UserData.headerString'])));
subpatternSeq=str2num(valueFromHeaderString('state.DMD.subpatternsString',...
  evalin('base',['AD0_e',num2str(epoch),'p2avg.UserData.headerString'])));

clearvars trialsMout

for i=1:length(traceNums)
    load([traceNums{i},'.mat'])
    trialsM(i,:) = eval([traceNums{i},'.data']);
end

clearvars AD0_* 

%% Filter out sweeps that are outside 25% range
clearvars sweepNum
if isempty(dir(sprintf('physParamsEpoch%d-passQC.mat',epoch)))
    [outliers,trialsMout] = filterSweepsQC(epoch,trialsM,10,0);
else
    load(sprintf('physParamsEpoch%d-passQC.mat',epoch),'sweepNum')
    load(sprintf('physParamsEpoch%d-failQC.mat',epoch),'outliers')
    trialsMout = trialsM(sweepNum,:);   
end

h=plot(trialsM','k');
hold on;plot(trialsMout');
title('gray traces are the sorted out sweeps')
set(gcf, 'Color', 'w');
savefig(gcf,'AfterSweepSortout');
close all

%% Pre-processing
% Calculate baseline current amplitude
base = mean(trialsMout(:,1:2999),2);
baseM = repmat(base,1,size(trialsMout,2));
trialsM_subtracted = trialsMout-baseM;
Fs = 10000;            % Sampling frequency    

% Low pass filter at 2kHz
LP=lowpass(trialsM_subtracted',2000,Fs);
% 
% % Notch filter
% d = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                'DesignMethod','butter','SampleRate',Fs);
% Notch = filtfilt(d,LP);   

% Smooth data using sgolay filter
yT=sgolayfilt(LP,5,27); % polynomial order of 5 and framelength of 27
y=yT';

% Median filter using 0.5ms window
y=movmedian(y,6,2);

% % If we don't pre-process
% y=trialsM_subtracted;
%% concatenate matrix for spot
clearvars tlist m v MAD med correctedSweeps
numPattern= length(patternSeq);

tlist = 300:100:(300+(numPattern-1)*100);
intv = [-30 70];  % time interval after stim to test
len = diff(intv)*10;

for ii = 1:length(tlist)
    clearvars rawSweeps preprocessedSweeps baseLoc offset
    t = tlist(ii); % in ms  
    baseLoc = mean(y(:,((t-30)*10+1):(t)*10),2);   % baseline calculated from 30ms before photostim.
  
    % remove outlier traces based on Raccess
    rawSweeps = trialsM(:,[(t+intv(1))*10+1:(t+intv(2))*10]);
    rawSweeps(outliers,:)=[];
    
    % baseline subtract again to make sure the outliers re the trace is centered
    % along zero
    preprocessedSweeps = y(:,[(t+intv(1))*10+1:(t+intv(2))*10])-repmat(baseLoc,1,len);
    % we do not need to remove outliers here again!! trialsMout already has
    % eliminated sweeps!!
    preprocessedSweeps_base = y(:,[(t+intv(1)-13)*10+1:(t+intv(2))*10])-repmat(baseLoc,1,len+130);
    
    
%     preprocessedSweeps = rawSweeps;   % if we use raw traces to do the
%     same analysis
    % baseline offset around stim artifact
    offset = mean(preprocessedSweeps(:,abs(intv(1))*10:(abs(intv(1))*10+50)),2);
    
%     % this is to capture spontaneous activity statistics
%     baseline = preprocessedSweeps(:,1:200)-offset;
%     
%     for base_i=1:size(baseline,1)
%         [Ebase(ii,base_i),Ebaseloc(ii,base_i)]=min(baseline(base_i,:));
%         
%         [Ibase(ii,base_i),Ibaseloc(ii,base_i)]=max(baseline(base_i,:));
%         
%     end
    
    % filter out sweeps where there was spontaenous event during photostim
    elim=find(abs(offset)>3);   % larger than 3pA is filtered out
    preprocessedSweeps(elim,:)=[];
    preprocessedSweeps_base(elim,:)=[];
    offset(elim)=[];
    
    % correct the raw traces with offset calculted over stim artifact
    correctedSweeps{ii} = preprocessedSweeps-offset;
    correctedSweeps_base{ii}= preprocessedSweeps_base-offset;
    m(ii,:) = mean(correctedSweeps{ii});
    v(ii,:) = var(correctedSweeps{ii},[],1);
    MAD(ii,:) = mad(correctedSweeps{ii},1);
    med(ii,:) = median(correctedSweeps{ii});
end

% % compare E and I size during baseline
% figure
% plot(-Ebase,Ibase,'ko')

%% plot all of them

figure
subplot(2,2,1)
x = linspace(intv(1),intv(2),length(m));
plot(x,m')
title('mean')
axis tight

subplot(2,2,2)
plot(x,v')
title('variance')
axis tight

subplot(2,2,3)
plot(x,med')
title('median')
axis tight

subplot(2,2,4)
plot(x,MAD')
title('median absoulte deviation')
axis tight

set(gcf,'color','w')
savefig(sprintf('hotspotFiltering_e%d',epoch))
close all
%% sort data into hotspots vs. non-hotspots
clearvars ipt pkWin_MAD
thres1=2;
% thres2=3; 
% hotspotsSUM=sum(MAD>thres1,2);
% hotspots=find(hotspotsSUM>50)';
hotspots=find(sum(isoutlier(MAD,'median',1),2)>50); % find spots that exceed 3 scaled MAD from the median for at least 5ms

%% Pick time window per hotspot based on change point analysis
clearvars pkWin ipt
for u=1:length(hotspots)
   clearvars timeRange 
    
%    timeRange=find(var(MAD)>(1.1*median(var(MAD))));   % threshold variance in MAD
   
%    timeRange= find(MAD(hotspots(u),:)>thres1); % simple threshold in MAD
%    timeRange= find(MAD(hotspots(u),:)>(median(MAD(hotspots(u),:)))); % simple threshold in MAD
   
   %change point analysis
%    ipt(u)=findchangepts(MAD(hotspots(u),:),'Statistic','std');
   ipt(u,:)=findchangepts([correctedSweeps{hotspots(u)}],'Statistic','rms','MaxNumChanges',2); % rms works better than std or mean
   ipt(u,1)=ipt(u,1)-cpaWind;   % widen the analysis window
   ipt(u,2)=ipt(u,2)+cpaWind;
%    if ~isempty(timeRange)
%         pkWin_MAD(u,:)=[min(timeRange),max(timeRange)];
%    end
   figure
   plot([correctedSweeps{hotspots(u)}]')
   hold on;
%    if ~isempty(timeRange)
%         vline(pkWin_MAD(u,1));vline(pkWin_MAD(u,2))
%    end
        vline(ipt(u,1),'k');vline(ipt(u,2),'k');
   set(gcf,'color','w')
   title(sprintf('spot #%d, cpa window increase =%.2fms',hotspots(u),cpaWind/10))
   savefig(sprintf('changePoint_spot%d',hotspots(u)))
   close all
end

%% Build a noise model for a given cell using parametric fitting

notHotspots = 1:size(MAD);
notHotspots(hotspots)=[];
notHotspotData=[];
for ui = 1:length(notHotspots)
    
    notHotspotData=[notHotspotData;correctedSweeps{notHotspots(ui)}];
    
end
% concatenate portion not containing stim artifact
notHotspotData=notHotspotData(:,355:end);
% [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(reshape(notHotspotData,1,[]));
notHotspotData=reshape(notHotspotData,1,[]);

figure;subplot(1,3,1);set(gcf,'color','w')
histogram(notHotspotData,'Normalization','pdf')
TF=isoutlier(notHotspotData);
hold on; histogram(notHotspotData(~TF),'Normalization','pdf');
pd1 = fitdist(notHotspotData(~TF)','Normal');
xrange=[min(notHotspotData),max(notHotspotData)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd1,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('not hotspots')

% let's build another noise model based on preStim data points 
preStimData=[];
for uj = 1:size(MAD)
    
    preStimData=[preStimData;correctedSweeps{uj}];
    
end
preStimData=preStimData(:,1:250);
preStimData=reshape(preStimData,1,[]);
subplot(1,3,2);histogram(preStimData,'Normalization','pdf')
TF2=isoutlier(preStimData);
hold on; histogram(preStimData(~TF2),'Normalization','pdf');
pd2 = fitdist(preStimData(~TF2)','Normal');
xrange=[min(preStimData),max(preStimData)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd2,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('pre-stim period')

subplot(1,3,3);
allNull=[notHotspotData(~TF),preStimData(~TF2)];
histogram(allNull,'Normalization','pdf')
pd3 = fitdist(allNull','Normal');
xrange=[min(allNull),max(allNull)];
hold on;plot(linspace(xrange(1),xrange(2)),pdf(pd3,linspace(xrange(1),xrange(2))),'LineWidth',2)
title('pre-stim period + not hotspot period')

%% save these parameters for a given cell
save('HotspotTimeWin','ipt','hotspots','cpaWind')
save('NoiseModel','pd1','pd2','pd3')
% there is always delay in MAD onset compared to actual jitter, hence I
% tried change point analysis in only the spots that crossed MAD threshold
close all

%% EPSC vs. IPSC peak comparison in the time window where change point analysis defined the onset and offset
clearvars avgEonly avgItoo pEPSC pEgivenI pI pIgivenE pboth

for u=1:length(hotspots)
    mkdir(sprintf('coreleaseAnalysis-spot%d',hotspots(u)))
    
    try
    clearvars E I Eloc Iloc adjInd_E Eavg adjInd_I Iavg tempMat tempMat_base E_null I_null baseline_u elim_u
    tempMat= correctedSweeps{hotspots(u)};
    tempMat_base=correctedSweeps_base{hotspots(u)};
    % find sweeps that cannot be used with stim onset as offset bc of
    % noisy signal or spontaneous even riding on top of it
%     baseline_u=mad(tempMat(:,1:(abs(intv(1))*10)),1,2);
    baseline_u = [mean(tempMat(:,1:100),2),mean(tempMat(:,(end-100):end),2)];
%     elim_u = union(find(sum(abs(baseline_u)>2,2)==1.5),find(diff(baseline_u,1,2)>5));
    elim_u = find(sum(abs(baseline_u)>1.5,2)==2);
    tempMat(elim_u,:)=[];
    tempMat_base(elim_u,:)=[];
    
    % calculate baseline peak sizes(-30ms to photostim onset)
    [E_null,E_nullLoc]=min(tempMat(:,1:(abs(intv(1))*10-5)),[],2);  % maybe should consider opening up this window to match the duration of point change interval 
    [I_null,I_nullLoc]=max(tempMat(:,1:(abs(intv(1))*10-5)),[],2);
%     figure;plot(-E_null,I_null,'o');axis([min(axis) max(axis) min(axis) max(axis)])
%     figure;plot(E_nullLoc,I_nullLoc,'o');axis([min(axis) max(axis) min(axis) max(axis)])
    % calculate EPSC and IPSC peak sizes within time window
    [E,Eloc]=min(tempMat(:,ipt(u,1):ipt(u,2)),[],2);
    [I,Iloc]=max(tempMat(:,ipt(u,1):ipt(u,2)),[],2);
    
    adjInd_E=Eloc+ipt(u,1)-1;
    adjInd_I=Iloc+ipt(u,1)-1;
    
    % avg around peak point
    for sweep_i = 1:size(tempMat,1)
        % to better estimate peak size post-stim period, average 0.5ms before and after peak and subtract from 10ms before 3ms of the peak
        base_beforepk(sweep_i) = mean(tempMat(sweep_i,(adjInd_E(sweep_i)-130):(adjInd_E(sweep_i)-30)),2);
        Eavg(sweep_i) = mean(tempMat(sweep_i,(adjInd_E(sweep_i)-5):(adjInd_E(sweep_i)+5)),2)-...
            base_beforepk(sweep_i);
        Iavg(sweep_i) = mean(tempMat(sweep_i,(adjInd_I(sweep_i)-5):(adjInd_I(sweep_i)+5)),2)-...
            base_beforepk(sweep_i);
        
        % to better estimate peak size pre-stim, baseline period
        base_beforepk_null(sweep_i) = mean(tempMat_base(sweep_i,(E_nullLoc(sweep_i)-130+130):(E_nullLoc(sweep_i)-30+130)),2);  % remember that tempMat_base has additional 130 points in the beginning
        E_nullavg(sweep_i) = mean(tempMat_base(sweep_i,(E_nullLoc(sweep_i)-5+130):(E_nullLoc(sweep_i)+5+130)),2)-...
            base_beforepk_null(sweep_i);
        I_nullavg(sweep_i) = mean(tempMat_base(sweep_i,(I_nullLoc(sweep_i)-5+130):(I_nullLoc(sweep_i)+5+130)),2)-...
            base_beforepk_null(sweep_i);

    end
    % range to show individual traces
    n0=1;
    n1=21;
    xx =linspace(intv(1),intv(2),length(tempMat));
    x0=0;
    y0=15;
    
    % Figure #1: example of peak detection across sweeps
    figure;set(gcf,'color','w')
    for (counter=n0:n1)
        plot(xx+x0*(counter-n0),tempMat(counter,:)+y0*(counter-n0));
        hold on;
        plot(adjInd_E(counter)/10+intv(1)+x0*(counter-n0),Eavg(counter)+y0*(counter-n0),'r*');
        plot(adjInd_I(counter)/10+intv(1)+x0*(counter-n0),Iavg(counter)+y0*(counter-n0),'b*');
        
        plot(xx+x0*(counter-n0),0*tempMat(counter,:)+y0*(counter-n0), 'color', 'black');
    end
    vline(ipt(u,1)/10+intv(1),':k')
    vline(ipt(u,2)/10+intv(1),':k')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig1',hotspots(u)));
    close all;
      
    % Figure #2: plot the EPSC vs. IPSC peak sizes along with null distribution 
    figure;set(gcf,'color','w')
%     plot(-E,I,'bo');hold on;
    % set threshold based on noise model - currently 2*sigma
    load('NoiseModel.mat','pd3')
    Ethres=-2*pd3.sigma;
    Ithres=2*pd3.sigma;
%     Ithres=1.96*sqrt(pd3.sigma);
%     % old version of setting threshold
%     Ethres=prctile(reshape(Ebase,1,[]),10); % set the threshold as lowest 10% percentile based on baseline fluctuations
%     Ithres=prctile(reshape(Ibase,1,[]),90);
    
    % calculate the probability of each quadrant
    pEPSC(u)=length(find(Eavg<Ethres))./length(Eavg);  % EPSC peak that is larger or equal to threshold 
    pEgivenI(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(find(Iavg>Ithres)); % conditional probability, given IPSC
    pEgivenNoI(u)=length(intersect(find(Eavg<Ethres),find(Iavg<=Ithres)))./length(find(Iavg<=Ithres)); % conditional probability, given no IPSC
    pIPSC(u)=length(find(Iavg>Ithres))./length(Iavg);  % EPSC peak that is larger or equal to threshold -- for a given data
    pIgivenE(u)=length(intersect(find(Iavg>Ithres),find(Eavg<Ethres)))./length(find(Eavg<Ethres)); % conditional probability, given EPSC
    pIgivenNoE(u)=length(intersect(find(Iavg>Ithres),find(Eavg>=Ethres)))./length(find(Eavg>=Ethres)); % conditional probability, given no EPSC
    pboth(u)=length(intersect(find(Eavg<Ethres),find(Iavg>Ithres)))./length(Eavg);
    % calculate correlation
    EICorr(u)=corr2(-Eavg,Iavg);
    EICorr_null(u)=corr2(-E_nullavg,I_nullavg);
    
    subplot(1,3,1)
    plot(-Eavg,Iavg,'ko')
    axis([min(axis) max(axis) min(axis) max(axis)])
    tempAx=axis();
    vline(-Ethres)
    hline(Ithres)
    xlabel('E peak (pA)')
    ylabel('I peak (pA)')
    
    text(max(axis)-4,2,sprintf('p(E)= %.2f',pEPSC(u)),'FontWeight','bold','Color','r')
    text(max(axis)-4,1,sprintf('p(E|I)= %.2f',pEgivenI(u)),'FontWeight','bold','Color','r')
    text(max(axis)-4,0,sprintf('p(E|no I)= %.2f',pEgivenNoI(u)),'FontWeight','bold','Color','r')
    text(0,max(axis)-4,sprintf('p(I)= %.2f',pIPSC(u)),'FontWeight','bold','Color','b')
    text(0,max(axis)-5,sprintf('p(I|E)= %.2f',pIgivenE(u)),'FontWeight','bold','Color','b')  
    text(0,max(axis)-6,sprintf('p(I|noE)= %.2f',pIgivenNoE(u)),'FontWeight','bold','Color','b')  
    text((max(axis)-min(axis))/1.7+min(axis),(max(axis)-min(axis))/1.7+min(axis),sprintf('p(E)p(I)= %.2f \n P(E|I)P(I)= %.2f \n P(I|E)P(E) = %.2f \n actual probability= %.2f',pEPSC(u).*pIPSC(u),pEgivenI(u).*pIPSC(u),pIgivenE(u).*pEPSC(u), pboth(u)),'FontWeight','bold')
    
    legend('avgpeak')
    title(sprintf('spot #%d',hotspots(u)))
    subplot(1,3,2)
    pEPSCn(u)=length(find(E_nullavg<Ethres))./length(E_nullavg);  % EPSC peak that is larger or equal to threshold
    pEgivenIn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg>Ithres)))./length(find(I_nullavg>Ithres)); % conditional probability, given IPSC
    pEgivenNoIn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg<=Ithres)))./length(find(I_nullavg<=Ithres)); % conditional probability, given no IPSC
    pIPSCn(u)=length(find(I_nullavg>Ithres))./length(I_nullavg);  % EPSC peak that is larger or equal to threshold
    pIgivenEn(u)=length(intersect(find(I_nullavg>Ithres),find(E_nullavg<Ethres)))./length(find(E_nullavg<Ethres)); % conditional probability, given EPSC
    pIgivenNoEn(u)=length(intersect(find(I_nullavg>Ithres),find(E_nullavg>=Ethres)))./length(find(E_nullavg>=Ethres)); % conditional probability, given no EPSC
    pbothn(u)=length(intersect(find(E_nullavg<Ethres),find(I_nullavg>Ithres)))./length(E_nullavg);
  
    sweep_idx=intersect(find(E_nullavg<Ethres),find(I_nullavg>abs(Ithres)));% changed from Ethres for both to Ithres
    plot(-E_nullavg,I_nullavg,'ko')
    % label the dots with sweep num
    for uu= 1:length(sweep_idx)
        text(-E_nullavg(sweep_idx(uu))+0.5,I_nullavg(sweep_idx(uu)),sprintf('%d',sweep_idx(uu)))
    end
    axis(tempAx)
    vline(-Ethres)
    hline(Ithres)
    xlabel('E peak (pA)')
    ylabel('I peak (pA)')
    text(max(axis)-4,2,sprintf('p(E)= %.2f',pEPSCn(u)),'FontWeight','bold','Color','r')
    text(max(axis)-4,1,sprintf('p(E|I)= %.2f',pEgivenIn(u)),'FontWeight','bold','Color','r')
    text(max(axis)-4,0,sprintf('p(E|no I)= %.2f',pEgivenNoIn(u)),'FontWeight','bold','Color','r')
    text(0,max(axis)-4,sprintf('p(I)= %.2f',pIPSCn(u)),'FontWeight','bold','Color','b')
    text(0,max(axis)-5,sprintf('p(I|E)= %.2f',pIgivenEn(u)),'FontWeight','bold','Color','b')   
    text(0,max(axis)-6,sprintf('p(I|noE)= %.2f',pIgivenNoEn(u)),'FontWeight','bold','Color','b')  
    text((max(axis)-min(axis))/1.7+min(axis),(max(axis)-min(axis))/1.7+min(axis),sprintf('p(E)p(I)= %.2f \n P(E|I)P(I)= %.2f \n P(I|E)P(E) = %.2f \n actual probability= %.2f',pEPSCn(u).*pIPSCn(u),pEgivenIn(u).*pIPSCn(u),pIgivenEn(u).*pEPSCn(u), pbothn(u)),'FontWeight','bold')
    
    legend('baseline peaks')
    title(sprintf('spot #%d',hotspots(u)))
    subplot(1,3,3)
    if ~isempty(sweep_idx)
        plot(xx,tempMat(sweep_idx,:)')
        hold on;
        plot(E_nullLoc(sweep_idx)/10+intv(1),E_nullavg(sweep_idx),'r*')
        plot(I_nullLoc(sweep_idx)/10+intv(1),I_nullavg(sweep_idx),'b*')
    end
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig2',hotspots(u)));
    close all;
     
    % Align by EPSC pk
    hasE=find(Eavg<Ethres)+n0-1;    % which sweep has E peak?, which adjustment of sweep counting shift
    alignRange=[-150,200];  % aligning 15ms before and 25ms after peak
    Eloc2=1+adjInd_E(hasE-n0+1);    % restore the shift
    subRange=tempMat((hasE-n0+1),:);    % restore the shift
    aligned=[];
    
    for i= 1:length(hasE)
        clearvars aligned_temp
            aligned_temp=subRange(i,(Eloc2(i)+alignRange(1)):(Eloc2(i)+alignRange(2)));
            aligned(i,:)=aligned_temp-mean(aligned_temp(1:-alignRange(1)));     % offset baseline
    end
    
    % Figure #3: probability densities
    pkData.sweepNum = 1:length(Eavg);
    pkData.Eavg= Eavg;
    pkData.Iavg = Iavg;
    pkData.hasE = Eavg<Ethres;
    pkData.hasI = Iavg>Ithres;
    pkData.E_nullavg=E_nullavg;
    pkData.I_nullavg=I_nullavg;
    pkData.Eloc = adjInd_E/10+intv(1);
    pkData.Iloc = adjInd_I/10+intv(1);
    
    figure; set(gcf,'color','w')
    clear fig3 
    fig3(1,1)=gramm('x',-pkData.Eavg,'color',pkData.hasI);
    fig3(1,1).stat_bin();
    fig3(1,1).set_title('EPSC amplitude')
    fig3(1,2)=gramm('x',-pkData.Eavg,'color',pkData.hasI);
    fig3(1,2).stat_density();
    fig3(1,2).set_title('EPSC amplitude')
    
    fig3(2,1)=gramm('x',pkData.Iavg,'color',pkData.hasE);
    fig3(2,1).stat_bin();
    fig3(2,1).set_title('IPSC amplitude')
    fig3(2,2)=gramm('x',pkData.Iavg,'color',pkData.hasE);
    fig3(2,2).stat_density();
    fig3(2,2).set_title('IPSC amplitude')
    fig3.draw();
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig3',hotspots(u)));
    close all;
    
    % Figure #4: E vs. I peaks across sweeps
    figure; set(gcf,'color','w')
    subplot(2,1,1);plot(Iavg,'bo')
    hold on; plot(Eavg,'ro');
    ylim([-50 50]);hline(-Ithres,'k');hline(Ithres,'k')
    legend('post-stim');
    title(sprintf('spot #%d',hotspots(u)));
    subplot(2,1,2);plot(I_nullavg,'bo');ylim([-50 50]);
    hold on; plot(E_nullavg,'ro')
    ylim([-50 50]);hline(-Ithres,'k');hline(Ithres,'k')
    legend('pre-stim')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig4',hotspots(u)));
    close all;
    
    % Figure #5: sweeps that contain EPSC pk, aligned to pk location
    figure;plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),aligned');set(gcf,'color','w')
    m_aligned = mean(aligned);
    v_aligned = var(aligned,[],1);
    MAD_aligned = mad(aligned,1);
    med_aligned = median(aligned);
    title('sweeps that contain EPSC peak, alignd to peak location')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig5',hotspots(u)));
    close all;
    
    % Figure #6: variance statistics of aligned traces that contain EPSC pk
    figure;set(gcf,'color','w')
    plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),[m_aligned;med_aligned;v_aligned;MAD_aligned]')
    legend('mean','median','variance','MAD')
    title('variance statistics of aligned traces that contain EPSC peak')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig6',hotspots(u)));
    close all;
    
    % check if hasE sweeps that also has IPSC
    hasItoo = intersect(find(Iavg>abs(Ithres))+n0-1,hasE);  % switch from Ethres for both to Ithres
    Eloc3=1+adjInd_E(hasItoo-n0+1);
    subRange3=tempMat((hasItoo-n0+1),:);    % restore the shift
    
    % Figure #7: plot the covariance of EPSC and IPSC pk timing ****
    clearvars idxIpreE
    figure;set(gcf,'color','w')
    Iloc2 = 1+adjInd_I(hasE-n0+1);    % restore the shift
    plot(Eloc2/10+intv(1),Iloc2/10+intv(1),'ko');
    hold on;
    plot(Eloc2/10+intv(1),(Iloc2-Eloc2)/10,'ro');
    axis([min(axis) max(axis) min(axis) max(axis)]) 
    line([min(axis) max(axis)],[min(axis) max(axis)],'color','r')
    idxIpreE=find(Iloc2<Eloc2); % find sweeps in which I precedes E
    % label the dots with sweep num
    for uv= 1:length(idxIpreE)
        text(Eloc2(idxIpreE(uv))/10+0.5,Iloc2(idxIpreE(uv))/10,sprintf('%d',hasE(idxIpreE(uv))-n0+1))
    end
    hline(0,'b:')
    xlabel('E timing(ms)')
    ylabel('I timing (ms)')
    legend(sprintf('peak location, sweeps that contain EPSC > 2pA, varE=%.1f, varI= %.1f',var(Eloc2/10),var(Iloc2/10)),...
        sprintf('difference of E-I, var=%.1f',var((Iloc2-Eloc2)/10)))
    title(sprintf('covariance of EPSC and IPSC peak timing in sweeps that contain E ,spot #%d',hotspots(u)))
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig7',hotspots(u)));
    close all;
    
    % Figure #8: plot the IPSC preceding EPSC sweeps with location of the
    % peaks
    figure;set(gcf,'color','w')
    ii=1;
    for (counter=(hasE(idxIpreE)-n0+1))     
        plot(xx+x0*(counter-n0),tempMat(counter,:)+y0*ii);
        hold on;
        plot(adjInd_E(counter)/10+intv(1)+x0*(counter-n0),Eavg(counter)+y0*ii+base_beforepk(counter),'r*');
        plot(adjInd_I(counter)/10+intv(1)+x0*(counter-n0),Iavg(counter)+y0*ii+base_beforepk(counter),'b*');
        plot(xx+x0*(counter-n0),0*tempMat(counter,:)+y0*ii, 'color', 'black');
        ii=ii+1;
    end
    vline(ipt(u,1)/10+intv(1),':k')
    vline(ipt(u,2)/10+intv(1),':k')
    title('IPSC preceding EPSC sweeps')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig8',hotspots(u)));
    close all;
    
    % these are the sweeps that have EPSC only
    hasEonly = setdiff(hasE,find(Iavg>abs(Ithres))+n0-1);   % switched from Ethres for both to Ithres
    Eloc4=1+adjInd_E(hasEonly-n0+1);
    subRange4=tempMat((hasEonly-n0+1),:);    % restore the shift

    aligned_Eonly=[];
    for iii= 1:length(hasEonly)
        clearvars aligned_Eonly_temp
        aligned_Eonly_temp=subRange4(iii,(Eloc4(iii)+alignRange(1)):(Eloc4(iii)+alignRange(2)));
        aligned_Eonly(iii,:)=aligned_Eonly_temp-mean(aligned_Eonly_temp(1:-alignRange(1))); 
    end
    
    % Figure #9: Separate traces into EPSC only vs. EPSC & IPSC with
    % alignment at EPSC pk
    figure;subplot(2,1,1);plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),aligned_Eonly');hold on; plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),mean(aligned_Eonly),'k','LineWidth',3)
    set(gcf,'color','w')
    title('EPSC only sweeps')
    aligned_EI=[];
    for iv= 1:length(hasItoo)
        clearvars aligned_EI_temp
        aligned_EI_temp=subRange3(iv,(Eloc3(iv)+alignRange(1)):(Eloc3(iv)+alignRange(2)));
        aligned_EI(iv,:)=aligned_EI_temp-mean(aligned_EI_temp(1:-alignRange(1))); 
    end
    subplot(2,1,2);plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),aligned_EI');hold on; plot(linspace(alignRange(1),alignRange(2),size(aligned,2)),mean(aligned_EI),'k','LineWidth',3)
    title('EPSC & IPSC containing sweeps')
    savefig(sprintf('coreleaseAnalysis-spot%d/Fig9',hotspots(u)));
    close all;
    
    % calculate the average E and I sizes
        avgEonly(u) = mean(Eavg(hasEonly-n0+1));
        pEonly(u) = length(hasEonly)./(length(hasEonly)+length(hasItoo));
        avgItoo(u) = mean(Eavg(hasItoo-n0+1));
        pItoo(u) = length(hasItoo)./(length(hasEonly)+length(hasItoo));
    % Plot the remainder of the plots only if there are enough EPSC only
    % sweeps! - for Figures 8 and 9    
    if size(aligned_Eonly,1)>3
        
        % Figure #10: EPSC amplitude quantil-quantile comparison of EPSC only sweeps vs both 
        figure;q1=qqplot(-Eavg(hasItoo-n0+1),-Eavg(hasEonly-n0+1));hold on; 

        q1(1).MarkerEdgeColor = [1 0 0];
        
        slope(u)=diff(q1(2).YData)./diff(q1(2).XData);
        intercept(u)=q1(2).YData(1)-slope(u).*q1(2).XData(1);
        gof(u)=goodnessOfFit(q1(1).YData,slope(u).*q1(1).XData+intercept(u),'MSE'); % calculate root mean squared error
        
        q2=qqplot(Iavg(hasItoo-n0+1),Iavg(hasEonly-n0+1));hold on; 
        set(gcf,'color','w')
        axis([min(axis) max(axis) min(axis) max(axis)])
        tempAx=axis();
        vline(-Ethres)
        hline(-Ethres)
        line([min(axis) max(axis)],[min(axis) max(axis)],'color','k')
        xlabel('Amplitude(pA)')
        ylabel('Amplitude (pA)')
        legend('Eavgpeak','IavgPeak')
        title('Q-Q plot of EPSC only vs both sweeps')
        savefig(sprintf('coreleaseAnalysis-spot%d/Fig10',hotspots(u)));
        close all;
        
        % do a permutation test to get a range
        lenX= length(hasItoo);
        lenY= length(hasEonly);
        totalN = lenX+lenY;
        H0data = [-Eavg(hasItoo-n0+1),-Eavg(hasEonly-n0+1)];

        nBoot=1000;
        rng default
        % Figure #11: plot permutation results for q-q plot
        figure; set(gcf,'color','w')
        clearvars n_X n_Y slope_perm intercept_perm gof_perm
        for k = 1:nBoot
            shuffledData = H0data(randperm(totalN));
            XStar = shuffledData(1:lenX);   
            YStar = shuffledData(lenX+1:end);
            n=qqplot(XStar,YStar);hold on;
            slope_perm(k)=diff(n(2).YData)./diff(n(2).XData);
            intercept_perm(k)=n(2).YData(1)-slope_perm(k).*n(2).XData(1);
            gof_perm(k)=goodnessOfFit(n(1).YData,slope_perm(k).*n(1).XData+intercept_perm(k),'MSE'); % calculate root mean squared error

            % save these variables later for plotting error as heatmap?
    %         n_X(k,:)=n(1).XData;
    %         n_Y(k,:)=n(1).YData;
        end
        myAlpha=0.05;
        idxLo = floor((myAlpha/2) * nBoot);    % index corresponding to lower bound
        idxHi = ceil((1-(myAlpha/2)) * nBoot); % index corresponding to upper bound
        sorted_slope = sort(slope_perm);
        sorted_intcpt = sort(intercept_perm);
        sorted_gof = sort(gof_perm);
        
        confInterval_slope = [sorted_slope(idxLo),sorted_slope(idxHi)];
        confInterval_intcpt = [sorted_intcpt(idxLo),sorted_intcpt(idxHi)];
        confInterval_gof = [sorted_gof(idxLo),sorted_gof(idxHi)];
        
        axis(tempAx)
        vline(-Ethres)
        hline(-Ethres)
        line([min(axis) max(axis)],[min(axis) max(axis)],'color','k')
        title('permutation test')
        savefig(sprintf('coreleaseAnalysis-spot%d/Fig11',hotspots(u)));
        close all;
        
        % Figure #12: showing distribution of parameters from permutation
        % test and comparing where the sample sits
        figure; set(gcf,'color','w')
        subplot(1,3,1)
        histogram(slope_perm);hold on; V1=vline(slope(u));V1.LineWidth=3;
        vline(confInterval_slope(1),'k');vline(confInterval_slope(2),'k')    
        title('slope distribution')
        
        subplot(1,3,2)
        histogram(intercept_perm);hold on; V2=vline(intercept(u));V2.LineWidth=3;
        vline(confInterval_intcpt(1),'k');vline(confInterval_intcpt(2),'k')
        title('intercept distribution')
        
        subplot(1,3,3)
        histogram(gof_perm);hold on; V3=vline(gof(u));V3.LineWidth=3;
        vline(confInterval_gof(1),'k');vline(confInterval_gof(2),'k')
        title('gof distribution')
        savefig(sprintf('coreleaseAnalysis-spot%d/Fig12',hotspots(u)));
        close all;
        
        clearvars Eloc2 Iloc2 Eloc3 Eloc4
    end
    catch error
        disp(sprintf('Hotspot #%d-failed',hotspots(u)))
    end
    
    % save all the meta-data
    save(sprintf('coreleaseAnalysis-spot%d/metadata',hotspots(u)),'pkData')

end

% Figure #13: avg Eonly and Itoo EPSC amplitudes 
figure;set(gcf,'color','w');
plot(-avgItoo,-avgEonly,'ko')
axis([0 50 0 50])
vline(-Ethres)
hline(-Ethres)
line([min(axis) max(axis)],[min(axis) max(axis)],'color','k')
xlabel('Amplitude (w IPSC)(pA)')
ylabel('Amplitude (w.o. IPSC)(pA)')
% label the dots with spot num
for uu= 1:length(hotspots)
    text(-avgItoo(uu)+1,-avgEonly(uu)+1,sprintf('%d',hotspots(uu)))
end
savefig('Fig13');
close all;

% Figure #13: probability of Eonly and Itoo for EPSC containing traces
figure;set(gcf,'color','w');
plot(pItoo,pEonly,'ko')
axis([0 1 0 1])
line([min(axis) max(axis)],[min(axis) max(axis)],'color','k')
xlabel('probability of EPSC with IPSC')
ylabel('probability of EPSC only')
% label the dots with spot num
for uu= 1:length(hotspots)
    text(pItoo(uu)+0.01,pEonly(uu)+0.01,sprintf('%d',hotspots(uu)))
end
savefig('Fig14');
close all;

% save these parameters for a given cell
save('statistics','pEPSCn','pEPSC','pEgivenIn','pEgivenI','pEgivenNoIn','pEgivenNoI','pIPSCn','pIPSC','pIgivenEn','pIgivenE','pIgivenNoEn','pIgivenNoE','EICorr','EICorr_null')


end % outer most loop