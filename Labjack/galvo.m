% Shun Li, 2023/08/13
% Modified from Paolo and Jay's code

% this script uses the DAQ with an input and an output
% input is coming from an arduino as analog input
% output will be used to turn ON an LED or a laser

%% Calibration of galvo and laser

clear; close all

% Ask user
answer = questdlg("Do you want to calibrate galvo & laser?",...
    "Calibration","Yes","No","Yes");
switch answer
    case "Yes"; calibration = 1;
    case "No"; calibration = 0;
end

if calibration == 1
    disp("----- Calibration started -----");
    d = daq('ni'); d.Rate = 20000;
    
    % Define input channel
    % opto_blue = addinput(d,'Galvo','ai0','Voltage');
    addinput(d,'Galvo','ai1','Voltage');
    
    % Define output channel
    % galvo_blue = addoutput(d,'Galvo','ao0','Voltage');
    addoutput(d,'Galvo','ao1','Voltage');
    
    write(d,0);
    
    waitdlg = warndlg('Press OK to finish calibration', 'Finish calibration');
    waitfor(waitdlg);
    write(d,2);
    disp("----- Calibration finished -----");
    clear; close all
end

%% Setup

clear; close all
addpath(addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Matlab')));

% Open a new daq session
d = daq('ni');
d.Rate = 20000;
d.DigitalTriggerTimeout = Inf; % time matlab waits for start trigger (sec)

% Define output channel
galvo_red = addoutput(d,'Galvo','ao1','Voltage');

% Defein trigger
addtrigger(d,"Digital","StartTrigger","External",'Galvo/PFI0');
d.NumDigitalTriggersPerRun = Inf;

% Create the pulse pattern
galvo_pattern = getGalvoPattern(pulse_width=5,pulse_freq=20,style="step");

%% Run session

waitdlg = warndlg('Press OK to start galvo control', 'Start session');
waitfor(waitdlg);

stop(d); flush(d); 
preload(d, galvo_pattern); 
start(d);

return

%% Function to generate galvo pattern
function galvo_pattern = getGalvoPattern(options)

    arguments
        options.pulse_width double = 10 % Time in ms
        options.pulse_freq double = 40
        options.stim_duration double = 500 % Time in ms
        options.samp_rate double = 20000
        options.off_voltage double = 2
        options.style string = "smooth"
    end

    if strcmp(options.style,"smooth")
        x_range = (1:options.samp_rate*options.stim_duration*0.001);
        wvf = (1 + sin((x_range)*options.pulse_freq*2*pi/options.samp_rate))*(options.off_voltage/2);
        galvo_pattern = wvf';
    elseif strcmp(options.style,"step")
        % Define pulse params
        pulse_width = 10 * 0.001 * options.samp_rate;
        pulse_freq = options.pulse_freq;
    
        nPulse = (options.stim_duration/1000) * pulse_freq;
        pulse_interval = ((1000/pulse_freq) - options.pulse_width) * 0.001 * options.samp_rate;
    
        if pulse_interval <= 0 && nPulse > 1
            pulse_interval = 10 * 0.001 * options.samp_rate;
            pulse_width = ((1000/pulse_freq)-10) * 0.001 * options.samp_rate;
            warning(strcat("Pulse interval is negative: reset pulse_interval=", ...
                num2str(pulse_interval/0.001/options.samp_rate), "ms and pulse_width=", ...
                num2str(pulse_width/0.001/options.samp_rate), "ms!"));
        end

        on_unit = zeros(pulse_width,1);
        off_unit = options.off_voltage * ones(pulse_interval,1);
    
        galvo_pattern = [];
        for i = 1:nPulse
            galvo_pattern = [galvo_pattern; on_unit; off_unit];
        end
    end

end
