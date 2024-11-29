function qc = getCellQC(headerString, options)

% Shun Li, 12/21/2023
% Modified from phUtil_rcCheck_Rin by Bernardo Sabatini
% Get QC (Rs, Rm, Cm) for each sweep

arguments
    headerString  char

    options.calculate logical = true % if false, just extract relevant value from headerString
    options.plot logical = false

    options.data double
    options.onset double
    options.pulseWidth double

    options.baselineWindow double
    options.Fs double = 10000
    options.currentClamp logical = false
    options.amplitude double = -5
    options.mode char {mustBeMember(options.mode,['peak', 'last10%', 'last30%', 'last50%'])} = 'last30%'

    options.rig string
end

%% Check input

if options.calculate && ~isfield(options,'data')
    error("Have to provide data if calculate from raw trace!");
end

%% Extract from headerString
if ~options.calculate
    qc.Rs = getHeaderValue(headerString,'state.phys.cellParams.rs0',convert=true);
    qc.Rm = getHeaderValue(headerString,'state.phys.cellParams.rm0',convert=true);
    qc.Cm = getHeaderValue(headerString,'state.phys.cellParams.cm0',convert=true);
    return
end

options.currentClamp = getHeaderValue(headerString,'state.phys.settings.currentClamp0');
options.Fs = getHeaderValue(headerString,'state.phys.settings.inputRate', convert=true); % ms

if isempty(getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck'))
    warning('No RC check detected!!!');
    qc.Rs = nan; qc.Rm = nan; qc.Cm = nan; 
    qc.Ibaseline = nan; qc.Ibaseline_std = nan; qc.Ibaseline_var = nan;
    qc.tau = nan; qc.Verror = nan;
    return
else
    options.amplitude = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='amplitude');
    options.onset = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='delay');
    options.pulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='pulseWidth');
    options.duration = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck', pulseVar='duration');
end

%% Adapted from Sanika's code

% Get important time points and I values for fitting decay 
pulsedt = 1000/options.Fs;
rcStartT = round((options.onset + 0.1)/pulsedt); % Takes about 0.1 ms for the voltage to be applied; need to avoid the capacititve component
rcEndT = round(rcStartT+(options.pulseWidth/pulsedt) - 1); % Stop before the next capacititve component
% 1 point before the end (previously 90% of the way as "-(options.pulseWidth/10)/pulsedt" instead of "- 1")

if options.amplitude > 0
    [Ipeak, IpeakT] = max(options.data(rcStartT:rcEndT)); % Save value of peak I and time index of peak I
else
    [Ipeak, IpeakT] = min(options.data(rcStartT:rcEndT));
end
IpeakT = IpeakT + rcStartT - 1;

% Get baseline current
if isfield(options,'baselineWindow')
    Ibaseline = median(options.data(options.baselineWindow)); % all baseline
    Ibaseline_std = std(options.data(options.baselineWindow));
    Ibaseline_var = var(options.data(options.baselineWindow));
else
    rcBaselineStart = round((rcStartT - 100)/pulsedt);
    options.baselineWindow = rcBaselineStart:rcStartT;
    Ibaseline = median(options.data(options.baselineWindow)); % baseline current before RC
    Ibaseline_std = std(options.data(options.baselineWindow));
    Ibaseline_var = var(options.data(options.baselineWindow));
end
qc.Ibaseline = Ibaseline;
qc.Ibaseline_std = Ibaseline_std; qc.Ibaseline_var = Ibaseline_var;

% Get steady state current
steadyStateT = (10/pulsedt); % In this case, the current has stabilized ~10 ms before the end (based on observation)
Isteadystate = median(options.data((rcEndT-steadyStateT):rcEndT)); % steadystate current during RC

% Taking 5 to 75% of the curve:
% Stop at 75% of the decay
Idecay75 = (Ipeak-Isteadystate)*.25 + Isteadystate; 
if options.amplitude > 0
    Idecay75T = find(options.data(IpeakT:rcEndT)<Idecay75, 1, 'first'); % T b/w Ipeak and Idecay75
else
    Idecay75T = find(options.data(IpeakT:rcEndT)>Idecay75, 1, 'first'); % T b/w Ipeak and Idecay75
end
% Starts after 5% of the decay 
Idecay05 = (Ipeak-Isteadystate)*.95 + Isteadystate; 
if options.amplitude > 0
    Idecay05T = find(options.data(IpeakT:rcEndT)<Idecay05, 1, 'first'); % T b/w Ipeak and Idecay05
else
    Idecay05T = find(options.data(IpeakT:rcEndT)>Idecay05, 1, 'first'); % T b/w Ipeak and Idecay05
end

%% Part 1: Get Rseries and Rmembrane  

decayStartT = IpeakT; % to enable back-extrapolation, or to include the first 5% of the decay
decayEndT = rcEndT; % to include steady state

% Alternatives, if desired:
% decayEndT = IpeakT + Idecay75T - 1; % Look at 0 to 75% of the decay
% decayStartT = IpeakT + Idecay05T; % Using the same decayEndT, look at 5% to 75% of the decay

% Perform the fit
decayTimeVectorTwoTerm = ((decayStartT:decayEndT) - decayStartT)*pulsedt;
decayDataVectorTwoTerm = options.data(decayStartT:decayEndT) - Isteadystate;
if options.amplitude < 0, decayDataVectorTwoTerm = -decayDataVectorTwoTerm; end

fitTypeTwoTerm = fittype('a*exp(-x/c) + b*exp(-x/d)', 'coeff', {'a', 'b', 'c', 'd'});
fitOptionsTwoTerm = fitoptions(fitTypeTwoTerm);
fitOptionsTwoTerm.Lower = [0 0 0 0];
fitOptionsTwoTerm.Upper = [abs(Ipeak-Isteadystate)*10, abs(Ipeak-Isteadystate)*10, options.duration, options.duration];
fitOptionsTwoTerm.StartPoint = [(Ipeak-Isteadystate)*0.8, (Ipeak-Isteadystate)*0.2, decayEndT*pulsedt/3, decayEndT*pulsedt/10];

fitFunctionTwoTerm = @(parameters, xData)parameters(1)*exp(-xData/parameters(3)) + parameters(2)*exp(-xData/parameters(4));
decayFitTwoTerm = fit(decayTimeVectorTwoTerm(:), decayDataVectorTwoTerm(:), fitTypeTwoTerm, fitOptionsTwoTerm);

%% Save Rseries and Rmembrane
% Find Rseries and Rmembrane

% final Rseries
Rseries = abs(options.amplitude)/(decayFitTwoTerm.a+decayFitTwoTerm.b) * 1000;
qc.Rs = Rseries;

% final Rmembrane
Rmembrane = (options.amplitude/(Isteadystate - Ibaseline) * 1000) - Rseries;
qc.Rm = Rmembrane;

qc.Rs_headerString = getHeaderValue(headerString,'state.phys.cellParams.rs0',convert=true);
qc.Rm_headerString = getHeaderValue(headerString,'state.phys.cellParams.rm0',convert=true);
qc.Cm_headerString = getHeaderValue(headerString,'state.phys.cellParams.cm0',convert=true);

%% Part 2: Get Tau with a single exponential fit

% Perform the fit
decayTimeVectorOneTerm = ((decayStartT:decayEndT) - decayStartT)*pulsedt;
decayDataVectorOneTerm = options.data(decayStartT:decayEndT) - Isteadystate;
if options.amplitude < 0, decayDataVectorOneTerm = -decayDataVectorOneTerm; end

fitTypeOneTerm = fittype('a*exp(-x/Tau)', 'coeff', {'a', 'Tau'});
fitOptionsOneTerm = fitoptions(fitTypeOneTerm);
fitOptionsOneTerm.Lower = [0 0];
fitOptionsOneTerm.Upper = [abs(Ipeak-Isteadystate)*10, options.duration];
fitOptionsOneTerm.StartPoint = [(Ipeak-Isteadystate)*0.8, rcEndT*pulsedt/3];

fitFunctionOneTerm = @(parameters, xData)parameters(1)*exp(-xData/parameters(2));
decayFitOneTerm = fit(decayTimeVectorOneTerm(:), decayDataVectorOneTerm(:), fitTypeOneTerm, fitOptionsOneTerm);
 
%% Save Rseries, Rmembrane, and Tau

% final Tau
Tau = decayFitOneTerm.Tau; % ms
qc.tau = Tau;

% final Cmembrane
Cmembrane = 1000 * Tau * ((1/Rseries) + (1/Rmembrane));
qc.Cm = Cmembrane;

% back up Rs, Rm, Cm (from one term fit)
Rseries_SingExp_noBackExtrap = abs(options.amplitude)/(decayFitOneTerm.a) * 1000; % why divide by decayFitOneTerm.a rather than Ipeak?
qc.Rs_SingleExp_noBackExtrap = Rseries_SingExp_noBackExtrap;

Rseries_SingExp = abs(options.amplitude/fitFunctionOneTerm([decayFitOneTerm.a, decayFitOneTerm.Tau], -.1))*1000;
qc.Rs_SingleExp = Rseries_SingExp;

Rmembrane_SingExp = (options.amplitude/(Isteadystate - Ibaseline) * 1000) - Rseries_SingExp;
qc.Rm_SingleExp = Rmembrane_SingExp;

Cmembrane_SingExp = 1000 * Tau * ((1/Rseries_SingExp) + (1/Rmembrane_SingExp));
qc.Cm_SingleExp = Cmembrane_SingExp;

% save the residuals -- must be over the same data range

%% Part 3: calculate voltage errors (added by Shun)

% V_error = Rs * I_baseline / 1000 (so its mV)
Verror = Rseries * Ibaseline / 1000;
qc.Verror = Verror;

%% Show plots

if options.plot
    initializeFig(0.8,0.8); tiledlayout(2,2); 
    
    nexttile; box off;
    plot(options.data, 'Color', [0 0.4470 0.7410]); hold on
    plot([rcStartT, IpeakT, decayStartT, decayEndT, rcEndT], [Ibaseline, Ipeak, Ipeak, Isteadystate, Isteadystate], '-o', 'Color', 'k');
    xlim([rcStartT-1000,rcEndT+1000]);
    
    nexttile; box off;
    plot(decayTimeVectorTwoTerm, decayDataVectorTwoTerm); hold on
    plot(decayTimeVectorTwoTerm, fitFunctionTwoTerm([decayFitTwoTerm.a, decayFitTwoTerm.b, decayFitTwoTerm.c, decayFitTwoTerm.d], decayTimeVectorTwoTerm));
    text(40, 700, "Two Term Fit");
    
    nexttile; box off;
    plot(decayTimeVectorOneTerm, decayDataVectorOneTerm); hold on
    plot(decayTimeVectorOneTerm, fitFunctionOneTerm([decayFitOneTerm.a, decayFitOneTerm.Tau], decayTimeVectorOneTerm));
    text(40, 700, "One Term Fit");
    
    nexttile; box off;
    plot(decayTimeVectorOneTerm, decayDataVectorOneTerm, 'Color', "k"); hold on
    plot(decayTimeVectorOneTerm, fitFunctionOneTerm([decayFitOneTerm.a, decayFitOneTerm.Tau], decayTimeVectorOneTerm), 'Color', "green");
    plot(decayTimeVectorTwoTerm, fitFunctionTwoTerm([decayFitTwoTerm.a, decayFitTwoTerm.b, decayFitTwoTerm.c, decayFitTwoTerm.d], decayTimeVectorTwoTerm), 'Color', "red");
    text(20, 700, "Data, (black), One Term (Green), Two Term (Red)");
end

end